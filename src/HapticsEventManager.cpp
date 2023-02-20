#include "HapticsEventManager.h"


using namespace midl;

float prev_stylus[3], curr_stylus[3], stylus[3], displacement[3], pivot[3];
  /******************************************************************************
  Factory method for creating an instance of the HapticsEventManager.
  ******************************************************************************/
HapticsEventManagerGeneric *HapticsEventManagerGeneric::Initialize()
{
	return new HapticsEventManager;
}

/******************************************************************************
Factory method for destroying an instance of the HapticsEventManager.
******************************************************************************/
void HapticsEventManagerGeneric::Delete(HapticsEventManagerGeneric *&pInterface)
{
	if (pInterface)
	{
		HapticsEventManager *pImp = static_cast<HapticsEventManager *>(pInterface);
		delete pImp;
		pInterface = 0;
	}
}

/******************************************************************************
HapticsEventManager Constructor.
******************************************************************************/
HapticsEventManager::HapticsEventManager() :
	m_hHD(HD_INVALID_HANDLE),
	m_hUpdateCallback(HD_INVALID_HANDLE),
	m_pHapticDeviceHT(0),
	m_pHapticDeviceGT(0),
	abcSketchManager(0),
	m_nCursorDisplayList(0)
{}

/*******************************************************************************
HapticsEventManager Destructor.
*******************************************************************************/
HapticsEventManager::~HapticsEventManager(){}

/*******************************************************************************
This is the main initialization needed for the haptic glue code.
*******************************************************************************/
void HapticsEventManager::Setup(ABCSketchManager *skManager)
{
	HDErrorInfo error;

	/* Intialize a device configuration. */
	m_hHD = hdInitDevice(HD_DEFAULT_DEVICE);
	if (HD_DEVICE_ERROR(error = hdGetError()))
	{
		std::cerr << error << std::endl;
		std::cerr << "Failed to initialize haptic device" << std::endl;
		std::cerr << "Press any key to quit." << std::endl;
		getchar();
		exit(-1);
	}

	/* Create the IHapticDevice instances for the haptic and graphic threads
	These interfaces are useful for handling the synchronization of
	state between the two main threads. */
	m_pHapticDeviceHT = IHapticDevice::create(
		IHapticDevice::HAPTIC_THREAD_INTERFACE, m_hHD);
	m_pHapticDeviceGT = IHapticDevice::create(
		IHapticDevice::GRAPHIC_THREAD_INTERFACE, m_hHD);

	/* Setup callbacks so we can be notified about events in the graphics
	thread. */
	m_pHapticDeviceGT->setCallback(
		IHapticDevice::MADE_CONTACT, madeContactCallbackGT, this);
	m_pHapticDeviceGT->setCallback(
		IHapticDevice::LOST_CONTACT, lostContactCallbackGT, this);
	m_pHapticDeviceGT->setCallback(
		IHapticDevice::BUTTON_1_UP, button1UpClickCallbackGT, this);
	m_pHapticDeviceGT->setCallback(
		IHapticDevice::BUTTON_1_DOWN, button1DownClickCallbackGT, this);
	m_pHapticDeviceGT->setCallback(
		IHapticDevice::BUTTON_2_UP, button2UpClickCallbackGT, this);
	m_pHapticDeviceGT->setCallback(
		IHapticDevice::BUTTON_2_DOWN, button2DownClickCallbackGT, this);
	m_pHapticDeviceGT->setCallback(
		IHapticDevice::DEVICE_ERROR, errorCallbackGT, this);

	hdEnable(HD_FORCE_OUTPUT);

	/* Schedule a haptic thread callback for updating the device every
	tick of the servo loop. */
	m_hUpdateCallback = hdScheduleAsynchronous(deviceUpdateCallback, this, HD_MAX_SCHEDULER_PRIORITY);

	/* Start the scheduler to get the haptic loop going. */
	hdStartScheduler();
	if (HD_DEVICE_ERROR(error = hdGetError()))
	{
		std::cerr << error << std::endl;
		std::cerr << "Failed to start scheduler" << std::endl;
		std::cerr << "Press any key to quit." << std::endl;
		getchar();
		exit(-1);
	}

	// GRAPHICS SETUP
	abcSketchManager =  skManager;
	diffuseShader.Initialize(".//Shaders//diffuseShader.vert", ".//Shaders//diffuseShader.frag");
	isForceOn = false;
}

/*******************************************************************************
Reverse the setup process by shutting down and destructing the services
used by the HapticsEventManager.
*******************************************************************************/
void HapticsEventManager::Cleanup()
{
	hdStopScheduler();

	if (m_hUpdateCallback != HD_INVALID_HANDLE)
	{
		hdUnschedule(m_hUpdateCallback);
		m_hUpdateCallback = HD_INVALID_HANDLE;
	}

	if (m_hHD != HD_INVALID_HANDLE)
	{
		hdDisableDevice(m_hHD);
		m_hHD = HD_INVALID_HANDLE;
	}

	IHapticDevice::destroy(m_pHapticDeviceGT);
	IHapticDevice::destroy(m_pHapticDeviceHT);

	glDeleteLists(m_nCursorDisplayList, 1);
}

/*******************************************************************************
This method will get called every tick of the graphics loop. It is primarily
responsible for synchronizing state with the haptics thread as well as
updating snap state.
*******************************************************************************/
void HapticsEventManager::UpdateState()
{
	/* Capture the latest state from the servoloop. */
	m_pHapticDeviceGT->beginUpdate(m_pHapticDeviceHT);
	m_pHapticDeviceGT->endUpdate(m_pHapticDeviceHT);
}

/*******************************************************************************
Uses the current OpenGL viewing transforms to determine a mapping from device
coordinates to world coordinates.
*******************************************************************************/
void HapticsEventManager::UpdateWorkspace()
{
	GLdouble modelview[16];
	GLdouble projection[16];
	GLint viewport[4];
		
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);
	
	/* Compute the transform for going from device coordinates to world
	coordinates based on the current viewing transforms. */
	hduMapWorkspaceModel(modelview, projection, m_workspaceXform);

	/* Compute the scale factor that can be applied to a unit sized object
	in world coordinates that will make it a particular size in pixels. */
	HDdouble screenTworkspace = hduScreenToWorkspaceScale(
		modelview, projection, viewport, m_workspaceXform);
	
	m_cursorScale = CURSOR_SIZE_PIXELS * screenTworkspace;

	hlMatrixMode(HL_TOUCHWORKSPACE);
	hluFitWorkspace(projection);
		
	/* Compute the updated camera position in world coordinates. */
	hduMatrix worldTeye(modelview);
	hduMatrix eyeTworld = worldTeye.getInverse();
	eyeTworld.multVecMatrix(hduVector3Dd(0, 0, 0), m_cameraPosWC);

	hdScheduleSynchronous(setDeviceTransformCallback, this, SYNCHRONIZE_STATE_PRIORITY);
}

/*******************************************************************************
Draws a 3D cursor using the current device transform and the workspace
to world transform.
*******************************************************************************/
void HapticsEventManager::RenderStylus()
{
	prev_stylus[0] = curr_stylus[0];
	prev_stylus[1] = curr_stylus[1];
	prev_stylus[2] = curr_stylus[2];
}


/******************************************************************************
Scheduler Callbacks
These callback routines get performed in the haptics thread.
******************************************************************************/

/******************************************************************************
This is the main haptic thread scheduler callback. It handles updating the
currently applied constraint.
******************************************************************************/
HDCallbackCode HDCALLBACK HapticsEventManager::deviceUpdateCallback(
	void *pUserData)
{	
	HapticsEventManager *pThis = static_cast<HapticsEventManager *>(pUserData);

	/* Force the haptic device to update its state. */
	pThis->m_pHapticDeviceHT->beginUpdate(0);

	IHapticDevice::IHapticDeviceState *pCurrentState =
		pThis->m_pHapticDeviceHT->getCurrentState();
	IHapticDevice::IHapticDeviceState *pLastState =
		pThis->m_pHapticDeviceHT->getLastState();

	/* Get the position of the device. */
	hduVector3Dd devicePositionLC = pCurrentState->getPosition();

	/* Get the gimbal angles of the device */
	HDdouble currentVelocity[3], currentAngularVelocity[3], jointAngle[3], gimbalAngle[3], currentPosition[3], transformMatrix[16];
	HDlong hapEncoder[6];

	hdGetDoublev(HD_CURRENT_VELOCITY, currentVelocity);
	hdGetDoublev(HD_CURRENT_ANGULAR_VELOCITY, currentAngularVelocity);
	hdGetDoublev(HD_CURRENT_JOINT_ANGLES, jointAngle);
	hdGetDoublev(HD_CURRENT_GIMBAL_ANGLES, gimbalAngle);
	hdGetDoublev(HD_CURRENT_POSITION, currentPosition);
	hdGetDoublev(HD_CURRENT_TRANSFORM, transformMatrix);
	hdGetLongv(HD_CURRENT_ENCODER_VALUES, hapEncoder);
	
	
	//float realstylus[] = { (float)devicePositionLC[0] , (float)devicePositionLC[1] , (float)devicePositionLC[2]  };

	float realVelocity[] = { (float)currentVelocity[0], (float)currentVelocity[1], (float)currentVelocity[2] };
	float realAngVelocity[] = { (float)currentAngularVelocity[0], (float)currentAngularVelocity[1], (float)currentAngularVelocity[2] };
	double realJointAngle[] = { jointAngle[0], jointAngle[1], jointAngle[2] };
	double realGimbalAngle[] = { gimbalAngle[0], gimbalAngle[1], gimbalAngle[2] };
	float realStylus[] = { (float)currentPosition[0] , (float)currentPosition[1] , (float)currentPosition[2]  };	
	float realMatrix[] = { (float)transformMatrix[0], (float)transformMatrix[1], (float)transformMatrix[2], (float)transformMatrix[3], (float)transformMatrix[4], (float)transformMatrix[5], (float)transformMatrix[6],
							(float)transformMatrix[7], (float)transformMatrix[8], (float)transformMatrix[9], (float)transformMatrix[10], (float)transformMatrix[11], (float)transformMatrix[12], (float)transformMatrix[13],
							(float)transformMatrix[14], (float)transformMatrix[15] };
	long realEncoder[6] = { hapEncoder[1], hapEncoder[2], hapEncoder[3], hapEncoder[4], hapEncoder[5], hapEncoder[6] };
		
	float forceVector[] = {0.0,0.0,0.0};

	/* Update the voxel model based on the new position data. */
	ABCSketchManager *skAPI = pThis->abcSketchManager;
	//Mesh *meshAPI = pThis->meshManager;
	
	skAPI->HapEncoderListen(realEncoder);
	skAPI->HapPositionListen(realStylus);
	skAPI->HapVelocityListen(realVelocity);
	skAPI->HapAngularVelocityListen(realAngVelocity);
	skAPI->HapJointAngleListen(realJointAngle);
	skAPI->HapGimbalAngleListen(realGimbalAngle);
	skAPI->HapMatrixListen(realMatrix);
	//skAPI->Listen(stylus);
	//bool isHapticsActive = skAPI->CheckHapticsForce();
	//skAPI->SetPivot(prev_stylus);
	//float pivot[3], 

	//float radius;

	/*displacement[0] = prev_stylus[0] - curr_stylus[0];
	displacement[1] = prev_stylus[1] - curr_stylus[1];
	displacement[2] = prev_stylus[2] - curr_stylus[2];
	float disp[3];
	Normalize3(displacement);*/

	//float rotaxis[3];
	//rotaxis[0] = displacement[1];
	//rotaxis[1] = -displacement[0];
	//rotaxis[2] = 0.0;

	//Normalize3(displacement);
	//float rotangle = -0.008*Norm3(displacement);
	//float vertex[3];	

	//float planePt[3]; skAPI->SetPlaneCenter(planePt);
	//float newOrigin[] = { 0.02,-0.31,0.25 }; //GL Coordinates
	//float planeNorm[] = { Zaxis[0], Zaxis[1], Zaxis[2] };

	//float planeDist = PointToPlaneDist(stylus, newOrigin, planeNorm);
	////cout << "planeDistance->" << planeDist << endl;
	//string option = skAPI->option;
	hduVector3Dd force;
	//bool isSketchOn = skAPI->isSketchOn;
	//bool isHapticsOn = skAPI->CheckHapticsForce();
	//cerr << "isHapticsOn" << isHapticsOn << endl;

	//if (isHapticsOn && option == "e" )
	//{
	//	float stiffness = 3.15;

	//	SubVectors3(curr_stylus, pivot, displacement);
	//	//Normalize3(displacement);
	//	
	//	float disp = Norm3(displacement);

	//	//cout << disp << endl;		
	//	
	//	forceVector[0] = 0.0; forceVector[1] = 0.0; forceVector[2] = stiffness * disp;
	//	force.set((HDdouble)forceVector[0], (HDdouble)forceVector[1], (HDdouble)forceVector[2]);
	//}

	//else if (isHapticsOn && option == "g")
	//{
	//	float radius[3], normal[3];
	//	float stiffness = 2.0;

	//	SubVectors3(stylus, newOrigin, radius);		

	//	//SubVectors3(pivot, Origin, radius);
	//	////float radii = PointToPlaneDist(pivot, Origin, Zaxis);
	//	////Normalize3(displacement);
	//	//float rad = Norm3(radius);

	//	//float radone = rad + 0.8;

	//	//float diff = normdist - rad;

	//	/*float temp = planeDist - rad;


	//	float radius2 = rad + 1.0;
	//	float radius1 = rad - 1.1;*/

	//	//if (planeDist >= radius1 && planeDist <= radius2)
	//	/*if (normdist <= radone)
	//	{*/
	//		
	//		//cerr << "within limits" << endl;
	//		//cerr << "Diff->" << diff << endl;
	//		//cerr << "temp->" << diff << endl;

	//		
	//	forceVector[0] = -stiffness *stylus[0]; forceVector[1] = -stiffness *stylus[1]; forceVector[2] = -stiffness *stylus[2];

	//	/*cerr << "OH: " << realstylus[0] << "  " << realstylus[1] << "  " << realstylus[2] << endl;
	//	cerr << "GL: " << stylus[0] << "  " << stylus[1] << "  " << stylus[2] << endl << endl;*/

	//	force.set((HDdouble)forceVector[0], (HDdouble)forceVector[1], (HDdouble)forceVector[2]);
	//}

	//else if (isHapticsOn && option == "l")
	//{
	//	float planeoffset[3], normal[3];
	//	float newOrigin[] = { 0.02,-0.31,0.25 }; //GL Coordinates
	//	float planeNorm[] = { Zaxis[0], Zaxis[1], Zaxis[2] };
	//	float planeDist = PointToPlaneDist(stylus, newOrigin, planeNorm);

	//	SubVectors3(pivot, newOrigin, planeoffset);
	//	float offset = PointToPlaneDist(pivot, newOrigin, planeNorm);
	//	float offsetMax = offset + 0.5;

	//	float diff = planeDist - offset;
	//	//cerr << "PD->" << planeDist << endl;
	//	
	//	float stiffness = 7.5;


	//	
	//	forceVector[0] = 0.0; forceVector[1] = 0.0; forceVector[2] = -stiffness*diff;
	//	force.set((HDdouble)forceVector[0], (HDdouble)forceVector[1], (HDdouble)forceVector[2]);
	//

	///*	else if (planeDist < dist1 && planeDist >= 0)
	//	{
	//		forceVector[0] = 0.0; forceVector[1] = 0.0; forceVector[2] = -stiffness * diff;
	//		force.set((HDdouble)forceVector[0], (HDdouble)forceVector[1], (HDdouble)forceVector[2]);
	//	}

	//	else if (planeDist < dist1 && planeDist < 0)
	//	{
	//		forceVector[0] = 0.0; forceVector[1] = 0.0; forceVector[2] = stiffness * diff;
	//		force.set((HDdouble)forceVector[0], (HDdouble)forceVector[1], (HDdouble)forceVector[2]);
	//	}*/
	//}
	////else if (isHapticsOn && option == "B" && planeDist < 0.18 && planeDist > -0.18  && stylus[1] > -2 && stylus[1] < 2 && stylus[0] > -1.5 && stylus[0] < 1.5 && isSketchOn)
	////{
	////	
	////	float stiffness = 9.5;
	////	float coeff = 0.41;
	////	forceVector[0] = coeff*disp[0]*stiffness*planeDist; forceVector[1] = coeff * disp[1] *stiffness*planeDist; forceVector[2] = -stiffness * planeDist;
	////	force.set((HDdouble)forceVector[0], (HDdouble)forceVector[1], (HDdouble)forceVector[2]);
	////	cerr << "Friction On" << endl;
	////	//cerr << "ForceX->" << (HDdouble)forceVector[0] << " " << "ForceY->" << (HDdouble)forceVector[1] << " " << "ForceZ->" << (HDdouble)forceVector[2] << endl;
	////}
	////else if (option == "U1" && planeDist < 0.06 && planeDist > -0.20 && stylus[1] > -2 && stylus[1] < 2 && stylus[0] > -1.5 && stylus[0] < 1.5 && isSketchOn)
	////{
	////	//meshAPI->GetAxisAngle(axis, -0.008*Norm3(displacement));
	////	float stiffness = 9.8;
	////	forceVector[0] = 0.0; forceVector[1] = 0.0; forceVector[2] = -stiffness*planeDist;
	////	force.set((HDdouble)forceVector[0], (HDdouble)forceVector[1], (HDdouble)forceVector[2]);
	////}
	////else if (isHapticsOn && option == "U1" && planeDist < 0.06 && planeDist > -0.20 && stylus[1] > -2 && stylus[1] < 2 && stylus[0] > -1.5 && stylus[0] < 1.5 && isSketchOn)
	////{
	////	//meshAPI->GetAxisAngle(axis, -0.008*Norm3(displacement));
	////	float stiffness = 9.5;
	////	float coeff = 0.41;
	////	forceVector[0] = coeff * disp[0] * stiffness*planeDist; forceVector[1] = coeff * disp[1] * stiffness*planeDist; forceVector[2] = -stiffness * planeDist;
	////	force.set((HDdouble)forceVector[0], (HDdouble)forceVector[1], (HDdouble)forceVector[2]);
	////}
	////else if (option == "U2" && planeDist < 0.20 && planeDist > 0.06 && stylus[1] > -2 && stylus[1] < 2 && stylus[0] > -1.5 && stylus[0] < 1.5 && isSketchOn)
	////{
	////	//meshAPI->GetAxisAngle(axis, -0.008*Norm3(displacement));
	////	float stiffness = 9.5;
	////	forceVector[0] = 0.0; forceVector[1] = 0.0; forceVector[2] = -stiffness*planeDist;
	////	force.set((HDdouble)forceVector[0], (HDdouble)forceVector[1], (HDdouble)forceVector[2]);
	////}
	////else if (isHapticsOn && option == "U2" && planeDist < 0.20 && planeDist > 0.06 && stylus[1] > -2 && stylus[1] < 2 && stylus[0] > -1.5 && stylus[0] < 1.5 && isSketchOn)
	////{
	////	//meshAPI->GetAxisAngle(axis, -0.008*Norm3(displacement));
	////	float stiffness = 9.5; 
	////	float coeff = 0.41;
	////	forceVector[0] = coeff * disp[0] * stiffness*planeDist; forceVector[1] = coeff * disp[1] * stiffness*planeDist; forceVector[2] = -stiffness * planeDist; forceVector[2] = -stiffness * planeDist;
	////	force.set((HDdouble)forceVector[0], (HDdouble)forceVector[1], (HDdouble)forceVector[2]);
	////}
	////else if (option == "C" && planeDist < 0.08 && planeDist > -0.08 && stylus[1] > -2 && stylus[1] < 2 && stylus[0] > -1.5 && stylus[0] < 1.5 && isSketchOn)
	////{
	////	//meshAPI->GetAxisAngle(axis, -0.008*Norm3(displacement));
	////	forceVector[0] = 0.0; forceVector[1] = 0.0; forceVector[2] = 1.5;
	////	force.set((HDdouble)forceVector[0], (HDdouble)forceVector[1], (HDdouble)forceVector[2]);
	////}
	////else if (isHapticsOn && option == "C" && planeDist < 0.08 && planeDist > -0.08 && stylus[1] > -2 && stylus[1] < 2 && stylus[0] > -1.5 && stylus[0] < 1.5 && isSketchOn)
	////{
	////	//meshAPI->GetAxisAngle(axis, -0.008*Norm3(displacement));
	////	float coeff = 0.81;
	////	forceVector[0] = coeff*1.5*disp[0]; forceVector[1] = coeff * 1.5*disp[1]; forceVector[2] = 1.5;
	////	force.set((HDdouble)forceVector[0], (HDdouble)forceVector[1], (HDdouble)forceVector[2]);
	////}
	//else 
	force.set((HDdouble)0.0, (HDdouble)0.0, (HDdouble)0.0);

	/*if (startSketch)
	force.set((HDdouble)forceVector[0], (HDdouble)forceVector[1], (HDdouble)forceVector[2]);
	else force.set((HDdouble)0.0, (HDdouble)0.0, (HDdouble)0.0);*/

	/*hduVector3Dd force;
	if (startSketch)
		force.set((HDdouble)forceVector[0], (HDdouble)forceVector[1], (HDdouble)forceVector[2]);
	else force.set((HDdouble)0.0, (HDdouble)0.0, (HDdouble)0.0);*/
	hdSetDoublev(HD_CURRENT_FORCE, force);

	//if (pCurrentState->isInContact())
	//{
	//	/* If currently in contact, use the freshest contact data. */
	//	pCurrentState->setContactData(pSnapAPI->getConstraint()->getUserData());
	//}
	//else if (pLastState->isInContact())
	//{
	//	/* If was in contact the last frame, use that contact data, since it
	//	will get reported to the event callbacks. */
	//	pCurrentState->setContactData(pLastState->getContactData());
	//}
	//else
	//{
	//	pCurrentState->setContactData((void *)-1);
	//}

	///* Transform result from world coordinates back to device coordinates. */
	//hduVector3Dd proxyPositionLC = pSnapAPI->getConstrainedProxy();
	//pCurrentState->setProxyPosition(proxyPositionLC);

	//hdGetDoublev(HD_CURRENT_TRANSFORM, pCurrentState->getProxyTransform());
	//pCurrentState->getProxyTransform()[12] = proxyPositionLC[0];
	//pCurrentState->getProxyTransform()[13] = proxyPositionLC[1];
	//pCurrentState->getProxyTransform()[14] = proxyPositionLC[2];

	//double kStiffness;
	//hdGetDoublev(HD_NOMINAL_MAX_STIFFNESS, &kStiffness);
	//kStiffness = hduMin(0.4, kStiffness);

	///* Compute spring force to attract device to constrained proxy. */
	//hduVector3Dd force = kStiffness * (proxyPositionLC - devicePositionLC);
	//hdSetDoublev(HD_CURRENT_FORCE, force);

	pThis->m_pHapticDeviceHT->endUpdate(0);

	return HD_CALLBACK_CONTINUE;
}

/******************************************************************************
Scheduler callback to set the workspace transform both for use in the graphics
thread and haptics thread.
******************************************************************************/
HDCallbackCode HDCALLBACK HapticsEventManager::setDeviceTransformCallback(
	void *pUserData)
{
	HapticsEventManager *pThis = static_cast<HapticsEventManager *>(pUserData);

	IHapticDevice::IHapticDeviceState *pStateGT =
		pThis->m_pHapticDeviceGT->getCurrentState();
	IHapticDevice::IHapticDeviceState *pStateHT =
		pThis->m_pHapticDeviceHT->getCurrentState();

	pStateGT->setParentCumulativeTransform(pThis->m_workspaceXform);
	pStateHT->setParentCumulativeTransform(pThis->m_workspaceXform);

	return HD_CALLBACK_DONE;
}

/******************************************************************************
Event Callbacks

These are event callbacks that are registered with the IHapticDevice
******************************************************************************/

/******************************************************************************
This handler gets called in the graphics thread whenever the device makes
contact with a constraint. Provide a visual cue (i.e. highlighting) to
accompany the haptic cue of being snapped to the point.
******************************************************************************/
void HapticsEventManager::madeContactCallbackGT(
	IHapticDevice::EventType event,
	const IHapticDevice::IHapticDeviceState * const pState,
	void *pUserData)
{
	HapticsEventManager *pThis = static_cast<HapticsEventManager *>(pUserData);

}

/******************************************************************************
This handler gets called in the graphics thread whenever the device loses
contact with a constraint. Provide a visual cue (i.e. highlighting) to
accompany the haptic cue of losing contact with the point.
******************************************************************************/
void HapticsEventManager::lostContactCallbackGT(
	IHapticDevice::EventType event,
	const IHapticDevice::IHapticDeviceState * const pState,
	void *pUserData)
{
	HapticsEventManager *pThis = static_cast<HapticsEventManager *>(pUserData);

}

/******************************************************************************
This handler gets called in the graphics thread whenever a button press is
detected. Interpret the click as a drilling on/off.
******************************************************************************/
void HapticsEventManager::button1UpClickCallbackGT(
	IHapticDevice::EventType event,
	const IHapticDevice::IHapticDeviceState * const pState,
	void *pUserData)
{
	HapticsEventManager *pThis = static_cast<HapticsEventManager *>(pUserData);

	ABCSketchManager *skAPI = pThis->abcSketchManager;
	/*skAPI->ManipulationOff();
	skAPI->RotationOff();
	skAPI->SketchOff();	*/
	skAPI->SaveLog();
}

void HapticsEventManager::button1DownClickCallbackGT(
	IHapticDevice::EventType event,
	const IHapticDevice::IHapticDeviceState * const pState,
	void *pUserData)
{
	HapticsEventManager *pThis = static_cast<HapticsEventManager *>(pUserData);

	//ABCSketchManager *skAPI = pThis->abcSketchManager;
	//skAPI->RotationOn();
	////pThis->RenderStylus();
	//pivot[0] = prev_stylus[0];
	//pivot[1] = prev_stylus[1];
	//pivot[2] = prev_stylus[2];

	//cout << pivot[0] << endl;
	//cout << pivot[1] << endl;
	//cout << pivot[2] << endl;	
}

void HapticsEventManager::button2UpClickCallbackGT(
	IHapticDevice::EventType event,
	const IHapticDevice::IHapticDeviceState * const pState,
	void *pUserData)
{
	HapticsEventManager *pThis = static_cast<HapticsEventManager *>(pUserData);

	ABCSketchManager *skAPI = pThis->abcSketchManager;
	//skAPI->ManipulationOff();
	
}

void HapticsEventManager::button2DownClickCallbackGT(
	IHapticDevice::EventType event,
	const IHapticDevice::IHapticDeviceState * const pState,
	void *pUserData)
{
	HapticsEventManager *pThis = static_cast<HapticsEventManager *>(pUserData);

	ABCSketchManager *skAPI = pThis->abcSketchManager;
	//skAPI->TranslationOn();
}


/******************************************************************************
This handler gets called to handle errors
******************************************************************************/
void HapticsEventManager::errorCallbackGT(
	IHapticDevice::EventType event,
	const IHapticDevice::IHapticDeviceState * const pState,
	void *pUserData)
{
	HapticsEventManager *pThis = static_cast<HapticsEventManager *>(pUserData);

	if (hduIsForceError(&pState->getLastError()))
	{
	}
	else
	{
		/* This is likely a more serious error, so just bail. */
		std::cerr << pState->getLastError() << std::endl;
		std::cerr << "Error during haptic rendering" << std::endl;
		std::cerr << "Press any key to quit." << std::endl;
		getchar();
		exit(-1);
	}
}

/******************************************************************************/

