#include <chrono>
#include <ctime>
#include "SketchManager.h"
#include "Renderer.h"


using namespace midl;
chrono::time_point<chrono::system_clock> timerON, timerOFF, currentTime;
__time64_t initTime, finalTime;

void ABCSketchManager::UpdateHapticEncoder(long *encoder)
{
	Tuple6l temp;

	temp.data[0] = encoder[0];
	temp.data[1] = encoder[1];
	temp.data[2] = encoder[2];
	temp.data[3] = encoder[3];
	temp.data[4] = encoder[4];
	temp.data[5] = encoder[5];

	rawEncoder.push_back(temp);
}

void ABCSketchManager::UpdateHapticStylus(float *stylus)
{
			Tuple3f temp;

			temp.data[0] = stylus[0];
			temp.data[1] = stylus[1];
			temp.data[2] = stylus[2];

			rawStylus.push_back(temp);
		//cerr << "Hap->"<< temp.data[0] << ", " << temp.data[1] << ", " << temp.data[2] << endl;
}

void ABCSketchManager::UpdateKinematicStylus(double *jointangle)
{
	/*gimbalPosition[0] = 0;
	gimbalPosition[1] = 0;
	gimbalPosition[2] = 0;*/

	float _gimbalPosition[3];

	_gimbalPosition[0] = 0.0;
	_gimbalPosition[1] = (-arm1 * cos(jointangle[1])) + (arm2*sin(jointangle[2]));
	_gimbalPosition[2] = (arm1 * sin(jointangle[1])) - (arm2*cos(jointangle[2]));

	Transformation T;
	T.SetYRotation(jointangle[0]);
	T.ApplyTo(_gimbalPosition);

	Tuple3f temp;

	temp.data[0] = _gimbalPosition[0];
	temp.data[1] = _gimbalPosition[1];
	temp.data[2] = _gimbalPosition[2];

	kinRawStylus.push_back(temp);

	//cerr << "Kinematic->" << temp.data[0] << ", " << temp.data[1] << ", " << temp.data[2] << endl;
}

void ABCSketchManager::UpdateHapticVelocity(float *velocity)
{
	Tuple3f temp;

	temp.data[0] = velocity[0];
	temp.data[1] = velocity[1];
	temp.data[2] = velocity[2];

	rawVelocity.push_back(temp);
}

void ABCSketchManager::UpdateHapticAngVelocity(float *angvelocity)
{
	Tuple3f temp;

	temp.data[0] = angvelocity[0];
	temp.data[1] = angvelocity[1];
	temp.data[2] = angvelocity[2];

	rawAngVelocity.push_back(temp);
}

void ABCSketchManager::UpdateHapticJointAngle(double *jointangle)
{
	Tuple3d temp;

	temp.data[0] = jointangle[0];
	temp.data[1] = jointangle[1];
	temp.data[2] = jointangle[2];

	rawJointAngle.push_back(temp);

	//cerr << "Hap->" << (180 / PI)*temp.data[1] << ", " << (180 / PI)*temp.data[2] << endl;
}

void ABCSketchManager::UpdateHapticGimbalAngle(double *gimbalangle)
{
	Tuple3d temp;

	temp.data[0] = gimbalangle[0];
	temp.data[1] = gimbalangle[1];
	temp.data[2] = gimbalangle[2];

	rawGimbalAngle.push_back(temp);

	

	/*cerr << "Kinematic->" << temp.data[0] << ", " << temp.data[1] << ", " << temp.data[2] << endl;*/
}

void ABCSketchManager::UpdateTransformMatrix(float *matrix)
{
	Tuple16f temp;

	temp.data[0] = matrix[0];
	temp.data[1] = matrix[1];
	temp.data[2] = matrix[2];
	temp.data[3] = matrix[3];
	temp.data[4] = matrix[4];
	temp.data[5] = matrix[5];
	temp.data[6] = matrix[6];
	temp.data[7] = matrix[7];
	temp.data[8] = matrix[8];
	temp.data[9] = matrix[9];
	temp.data[10] = matrix[10];
	temp.data[11] = matrix[11];
	temp.data[12] = matrix[12];
	temp.data[13] = matrix[13];
	temp.data[14] = matrix[14];
	temp.data[15] = matrix[15];	

	rawTransform.push_back(temp);
}

void ABCSketchManager::SetDrillingStatus(int status)
{
	if (status == -1)
	{
		isPreparing = true; //This is true
		isApproaching = false;
		isPullBack = false;
		isDrillingCortex1 = false;
		isDrillingCortex2 = false;

		cerr << "\n Preparing the drill...." << endl;
	}
	else if (status == 0)
	{
		isPreparing = false;
		isApproaching = true; //This is True
		isPullBack = false;
		isDrillingCortex1 = false;
		isDrillingCortex2 = false;

		cerr << "\n Approaching the bone..." << endl;
	}
	else if (status == 1)
	{
		isPreparing = false;
		isApproaching = false;
		isPullBack = false;
		isDrillingCortex1 = true; //This is true
		isDrillingCortex2 = false;

		cerr << "\n Drilling Through Cortex 1..." << endl;
	}
	else if (status == 2)
	{
		isPreparing = false;
		isApproaching = false;
		isPullBack = false;
		isDrillingCortex1 = false;
		isDrillingCortex2 = true; //This is true

		cerr << "\n Drilling Through Cortex 2...." << endl;
	}
	else if (status == 3)
	{
		isPreparing = false;
		isApproaching = false;
		isPullBack = true; //This is true
		isDrillingCortex1 = false;
		isDrillingCortex2 = false; 

		cerr << "\n Pulling Back the Drill...." << endl;
	}
	else if (status == 4)
	{
		isPullBack = false; //This is true
		isDrillingCortex1 = false;
		isDrillingCortex2 = false;

		cerr << "\n Drilling is Completed...." << endl;
	}
}

void ABCSketchManager::UpdateDrillingStatus()
{
	if (isPreparing) drillStatus.push_back(-1);
	else if (isApproaching) drillStatus.push_back(0);
	else if (isDrillingCortex1) drillStatus.push_back(1);
	else if (isDrillingCortex2) drillStatus.push_back(2);
	else if (isPullBack) drillStatus.push_back(3);
	else if (!isPullBack && !isDrillingCortex1) drillStatus.push_back(4);
}

ABCSketchManager::ABCSketchManager()
{
	isRecording = false;
	isPreparing = false;
	isApproaching = false;
	isPullBack = false; 
	isDrillingCortex1 = false;
	isDrillingCortex2 = false;
	isReset = true;
}

ABCSketchManager::~ABCSketchManager()
{
	rawEncoder.clear();
	rawStylus.clear();
	rawVelocity.clear();
	rawAngVelocity.clear();
	rawJointAngle.clear();
	rawGimbalAngle.clear();
	rawTransform.clear();

	isRecording = false;
	isPreparing = false;
	isApproaching = false;
	isPullBack = false;
	isDrillingCortex1 = false;
	isDrillingCortex2 = false;
	isReset = true;
}

void ABCSketchManager::InitShaders()
{
	diffuseShader.Initialize(".//Shaders//diffuseShader.vert", ".//Shaders//diffuseShader.frag");
	textureShader.Initialize(".//Shaders//colorShader.vert", ".//Shaders//colorShader.frag");
}

void ABCSketchManager::Update()
{
	if (isRecording)
	{
		currentTime = chrono::system_clock::now();
		chrono::duration<float> _timeElapsed = currentTime - timerON;
		timeElapsed.push_back(_timeElapsed.count());

		UpdateHapticEncoder(hapcurr_encoder);
		UpdateHapticStylus(hapcurr_stylus);
		UpdateKinematicStylus(hapcurr_joiangle);
		UpdateHapticVelocity(hapcurr_velocity);
		UpdateHapticAngVelocity(hapcurr_angvelocity);
		UpdateHapticJointAngle(hapcurr_joiangle);
		UpdateHapticGimbalAngle(hapcurr_gimangle);
		UpdateTransformMatrix(hap_matrix);
		UpdateDrillingStatus();
	}
}

void ABCSketchManager::HapEncoderListen(long *encoder)
{
	hapcurr_encoder[0] = encoder[0];
	hapcurr_encoder[1] = encoder[1];
	hapcurr_encoder[2] = encoder[2];
}

void ABCSketchManager::HapPositionListen(float *stylus)
{
	hapcurr_stylus[0] = stylus[0];
	hapcurr_stylus[1] = stylus[1];
	hapcurr_stylus[2] = stylus[2];
}

void ABCSketchManager::HapVelocityListen(float *velocity)
{
	hapcurr_velocity[0] = velocity[0];
	hapcurr_velocity[1] = velocity[1];
	hapcurr_velocity[2] = velocity[2];
}

void ABCSketchManager::HapAngularVelocityListen(float *angvelocity)
{
	hapcurr_angvelocity[0] = angvelocity[0];
	hapcurr_angvelocity[1] = angvelocity[1];
	hapcurr_angvelocity[2] = angvelocity[2];

	//cerr << angvelocity[0] << ", " << angvelocity[1] << ", " << angvelocity[2] << endl;
}

void ABCSketchManager::HapJointAngleListen(double *jointangle)
{
	hapcurr_joiangle[0] = jointangle[0];
	hapcurr_joiangle[1] = jointangle[1];
	hapcurr_joiangle[2] = jointangle[2];
}


void ABCSketchManager::HapGimbalAngleListen(double *gimbalangle)
{
	hapcurr_gimangle[0] = gimbalangle[0];
	hapcurr_gimangle[1] = gimbalangle[1];
	hapcurr_gimangle[2] = gimbalangle[2];
}

void ABCSketchManager::HapMatrixListen(float *matrix)
{
	hap_matrix[0] = matrix[0];
	hap_matrix[1] = matrix[1];
	hap_matrix[2] = matrix[2];
	hap_matrix[3] = matrix[3];
	hap_matrix[4] = matrix[4];
	hap_matrix[5] = matrix[5];
	hap_matrix[6] = matrix[6];
	hap_matrix[7] = matrix[7];
	hap_matrix[8] = matrix[8];
	hap_matrix[9] = matrix[9];
	hap_matrix[10] = matrix[10];
	hap_matrix[11] = matrix[11];
	hap_matrix[12] = matrix[12];
	hap_matrix[13] = matrix[13];
	hap_matrix[14] = matrix[14];
	hap_matrix[15] = matrix[15];
}

void ABCSketchManager::Reset()
{
	isRecording = false;
	isPreparing = false;
	isApproaching = false;
	isPullBack = false;
	isDrillingCortex1 = false;
	isDrillingCortex2 = false;
	//userID.clear();
	boneType.clear();
	//trial.clear();

	/*cerr << "\n Please enter the User ID ?" << endl;
	cin >> userID;
*/
	cerr << "\n Please enter the bone variant ? (OB/YB)" << endl;
	cin >> boneType;

	//cerr << "\n Please enter the trial number ?" << endl;
	//cin >> trial;
}


void ABCSketchManager::StartTimer()
{
	timerON = chrono::system_clock::now();
	initTime = chrono::system_clock::to_time_t(timerON);

	cerr << "\n Recording Started at->" << ctime(&initTime) << endl;
}

void ABCSketchManager::RecordData()
{
	StartTimer();
	isRecording = true;
	isReset = false;
	//cerr << "Recording Started" << endl;
}


void ABCSketchManager::SaveLog()
{
	isRecording = false;	
	//SaveKinematicLog();

	ofstream myfile;

	if (boneType == "OB" || boneType == "ob")
	{
		OBtrial++;
		myfile.open("PO_User" + to_string(userID) + "_" + boneType + "_" + to_string(OBtrial) + ".txt");
	}
	else if (boneType == "YB" || boneType == "yb")
	{
		YBtrial++;
		myfile.open("PO_User" + to_string(userID) + "_" + boneType + "_" + to_string(YBtrial) + ".txt");
	}
	

	for (int i = 0; i < rawStylus.size(); i++)
	{
		myfile << drillStatus[i] <<" "<< timeElapsed[i]<<" "<<rawEncoder[i].data[0] << " " << rawEncoder[i].data[1] << " " << rawEncoder[i].data[2] << " " << rawEncoder[i].data[3] << " " << rawEncoder[i].data[4] << " " << rawEncoder[i].data[5] << " "
			<< rawStylus[i].data[0] << " " << rawStylus[i].data[1] << " " << rawStylus[i].data[2] << " "
			<< rawVelocity[i].data[0] << " " << rawVelocity[i].data[1] << " " << rawVelocity[i].data[2] << " "
			<< rawAngVelocity[i].data[0] << " " << rawAngVelocity[i].data[1] << " " << rawAngVelocity[i].data[2] << " "
			<< rawJointAngle[i].data[0] << " " << rawJointAngle[i].data[1] << " " << rawJointAngle[i].data[2] << " "
			<< rawGimbalAngle[i].data[0] << " " << rawGimbalAngle[i].data[1] << " " << rawGimbalAngle[i].data[2] << " "
			<< rawTransform[i].data[0] << " " << rawTransform[i].data[1] << " " << rawTransform[i].data[2] << " " << rawTransform[i].data[3] << " " << rawTransform[i].data[4] << " " << rawTransform[i].data[5] << " "
			<< rawTransform[i].data[6] << " " << rawTransform[i].data[7] << " " << rawTransform[i].data[8] << " " << rawTransform[i].data[9] << " " << rawTransform[i].data[10] << " " << rawTransform[i].data[11] << " "
			<< rawTransform[i].data[12] << " " << rawTransform[i].data[13] << " " << rawTransform[i].data[14] << " " << rawTransform[i].data[15] << endl;
	}
	SaveTimeStamp();

	//cerr << "Data Logged Successfully" << endl;



	rawEncoder.clear();
	rawStylus.clear();
	rawVelocity.clear();
	rawAngVelocity.clear();
	rawJointAngle.clear();
	rawGimbalAngle.clear();
	rawTransform.clear();
	timeElapsed.clear();
	drillStatus.clear();

	myfile.close();
}

//void ABCSketchManager::SaveKinematicLog()
//{
//	ofstream myfile;
//	myfile.open("KinLog_" + to_string(userID) + "_" + to_string(trial) + ".txt");
//
//	for (int i = 0; i < rawStylus.size(); i++)
//	{
//		myfile << timeElapsed[i] << kinRawStylus[i].data[0] << " " << kinRawStylus[i].data[1] << " " << kinRawStylus[i].data[2] << " "
//			<< rawJointAngle[i].data[0] << " " << rawJointAngle[i].data[1] << " " << rawJointAngle[i].data[2] << endl;
//	}
//
//	timeElapsed.clear();
//	kinRawStylus.clear();
//	rawJointAngle.clear();
//}

void ABCSketchManager::SaveTimeStamp()
{
	timerOFF = chrono::system_clock::now();

	
	finalTime = chrono::system_clock::to_time_t(timerOFF);

	chrono::duration<float> timeTAKEN = timerOFF - timerON;

	ofstream myfile;
	if (boneType == "OB" || boneType == "ob")
	{
		myfile.open("PT_User" + to_string(userID) + "_" + boneType + "_" + to_string(OBtrial) + ".txt");
	}
	else if (boneType == "YB" || boneType == "yb")
	{
		myfile.open("PT_User" + to_string(userID) + "_" + boneType + "_" + to_string(YBtrial) + ".txt");
	}

	myfile << ctime(&initTime) << endl; 
	myfile << ctime(&finalTime) << endl;
	//myfile << timeTAKEN.count() << endl;

	cerr << "\n Recording Stopped at->" << ctime(&finalTime) << endl;
	cerr << "\n The task took->" << timeTAKEN.count() << " " << "seconds" << endl;

	myfile.close();

	Reset();
	//cerr << "\n Hit Esc to close OR Hit X/x to reset" << endl;
}
