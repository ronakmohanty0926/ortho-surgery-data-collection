#pragma once

#include "Core.h"
#include "HapticsEventManager.h"

namespace midl {

	class ABCSketchManager
	{
	private:
		//Shaders
		Shader diffuseShader, textureShader;

		//Haptic Device Variables
		float hapcurr_stylus[3], hapcurr_velocity[3], hapcurr_angvelocity[3], hap_matrix[16];
		long hapcurr_encoder[6];
		double hapcurr_gimangle[3], hapcurr_joiangle[3];
		
		//Haptic Device Physical Measurements
		float arm1, arm2, gimbalPosition[3];

		//Clock Variables
		int _timerValue;
		bool timerSet;

		bool isPreparing;
		bool isApproaching;
		bool isDrillingCortex1;
		bool isDrillingCortex2;
		bool isPullBack;

		//Haptic Device Data Storage Variables
		vector<float> timeElapsed;
		vector<int> drillStatus;
		vector<Tuple6l> rawEncoder;
		vector<Tuple3f> rawStylus;
		vector<Tuple3f> kinRawStylus;
		vector<Tuple3f> rawVelocity;
		vector<Tuple3f> rawAngVelocity;
		vector<Tuple3d> rawJointAngle;
		vector<Tuple3d> rawGimbalAngle;
		vector<Tuple16f> rawTransform;

		//Calibration Data Storage Variables
		//vector<Tuple3f> gimbalPosition;
		
		//Haptic Device Data Update Variables
		void UpdateHapticEncoder(long *encoder);
		void UpdateHapticStylus(float *stylus);
		void UpdateKinematicStylus(double *jointangle);
		void UpdateHapticVelocity(float *velocity);
		void UpdateHapticAngVelocity(float *angvelocity);
		void UpdateHapticJointAngle(double *jointangle);
		void UpdateHapticGimbalAngle(double *gimbalangle);
		void UpdateTransformMatrix(float *matrix);
		void UpdateDrillingStatus();

		void StartTimer();
		
		void SaveTimeStamp();
		void SaveKinematicLog();
	
	public:		
		ABCSketchManager();
		~ABCSketchManager();

		//string username;
		string boneType;
		int OBtrial =  -1;
		int YBtrial = -1;
		int userID;
		bool isRecording;
		bool isReset;

		//Shader Function
		void InitShaders();		
	
		//Global Update Functions
		void Update();
		void HapEncoderListen(long *encoder);
		void HapPositionListen(float *stylus);
		void HapVelocityListen(float *velocity);
		void HapAngularVelocityListen(float *angvelocity);
		void HapJointAngleListen(double *jointangle);
		void HapGimbalAngleListen(double *gimbalangle);
		void HapMatrixListen(float *matrix);
		void Reset();

		//Data Saving
		void RecordData();
		void SaveLog();	
		void SetDrillingStatus(int status);

	};
}