// Copyright © 2008-2013 Pioneer Developers. See AUTHORS.txt for details
// Licensed under the terms of the GPL v3. See licenses/GPL-3.txt

#include "libs.h"
#include "Pi.h"
#include "OculusRift.h"

#include "LibOVR/Include/OVR.h"
#include "LibOVR/Include/OVRVersion.h"
#include "LibOVR/Src/Util/Util_Render_Stereo.h"
using namespace OVR;
using namespace OVR::Util::Render;

// OculusRiftImplemetation:
// 
// Concrete implementation of the Oculus Rift / LibOVR management - non-virtual.
// Purpose is to isolate the LibOVR from Pioneer and allow for changes to the API & SDK without affecting Pioneer.
class OculusRiftImplemetation : public MessageHandler {
public:
	OculusRiftImplemetation() : pManager(nullptr), pSensor(nullptr), pHMD(nullptr)
	{
		Width  = 1280;
		Height = 800;

		// Initializes LibOVR. This LogMask_All enables maximum logging.
		// Custom allocator can also be specified here.
		System::Init(Log::ConfigureDefaultLog(LogMask_All));
		if(!System::IsInitialized())
			return;

		// *** Oculus HMD & Sensor Initialization

		// Create DeviceManager and first available HMDDevice from it.
		// Sensor object is created from the HMD, to ensure that it is on the
		// correct device.

		pManager = DeviceManager::Create();
		if(!pManager)
			return;

		// We'll handle it's messages in this case.
		pManager->SetMessageHandler(this);


		int         detectionResult = IDCONTINUE;
		const char* detectionMessage;

		do 
		{
			// Release Sensor/HMD in case this is a retry.
			if(pSensor) { delete pSensor; pSensor=nullptr; }
			if(pHMD) { delete pHMD; pHMD=nullptr; }

			pHMD = pManager->EnumerateDevices<HMDDevice>().CreateDevice();
			if (nullptr!=pHMD)
			{
				pSensor = pHMD->GetSensor();

				// This will initialize HMDInfo with information about configured IPD,
				// screen size and other variables needed for correct projection.
				// We pass HMD DisplayDeviceName into the renderer to select the
				// correct monitor in full-screen mode.
				if (pHMD->GetDeviceInfo(&HMDInfo))
				{            
					SConfig.SetHMDInfo(HMDInfo);
				}
			}
			else
			{            
				// If we didn't detect an HMD, try to create the sensor directly.
				// This is useful for debugging sensor interaction; it is not needed in
				// a shipping app.
				pSensor = pManager->EnumerateDevices<SensorDevice>().CreateDevice();
			}


			// If there was a problem detecting the Rift, display appropriate message.
			detectionResult  = IDCONTINUE;        

			if (!pHMD && !pSensor)
				detectionMessage = "Oculus Rift not detected.";
			else if (!pHMD)
				detectionMessage = "Oculus Sensor detected; HMD Display not detected.";
			else if (!pSensor)
				detectionMessage = "Oculus HMD Display detected; Sensor not detected.";
			else if (HMDInfo.DisplayDeviceName[0] == '\0')
				detectionMessage = "Oculus Sensor detected; HMD display EDID not detected.";
			else
				detectionMessage = nullptr;

			if (detectionMessage)
			{
				/*String messageText(detectionMessage);
				messageText += "\n\n"
							   "Press 'Try Again' to run retry detection.\n"
							   "Press 'Continue' to run anyway.";
							   */
				detectionResult = IDCONTINUE;//::MessageBoxA(0, messageText.ToCStr(), "Oculus Rift Detection", MB_CANCELTRYCONTINUE|MB_ICONWARNING);

				//if (detectionResult == IDCANCEL)
				//	return;
			}

		} while (detectionResult != IDCONTINUE);

    
		if (HMDInfo.HResolution > 0)
		{
			Width  = HMDInfo.HResolution;
			Height = HMDInfo.VResolution;
		}


		//if (!setupWindow())
		//	return 1;
    
		if (pSensor)
		{
			// We need to attach sensor to SensorFusion object for it to receive 
			// body frame messages and update orientation. SFusion.GetOrientation() 
			// is used in OnIdle() to orient the view.
			SFusion.AttachToSensor(pSensor);
			SFusion.SetDelegateMessageHandler(this);
			SFusion.SetPredictionEnabled(true);
		}

    	// *** Configure Stereo settings.

		SConfig.SetFullViewport(Viewport(0,0, Width, Height));
		SConfig.SetStereoMode(Stereo_LeftRight_Multipass);

		// Configure proper Distortion Fit.
		// For 7" screen, fit to touch left side of the view, leaving a bit of invisible
		// screen on the top (saves on rendering cost).
		// For smaller screens (5.5"), fit to the top.
		if (HMDInfo.HScreenSize > 0.0f)
		{
			if (HMDInfo.HScreenSize > 0.140f) // 7"
				SConfig.SetDistortionFitPointVP(-1.0f, 0.0f);
			else
				SConfig.SetDistortionFitPointVP(0.0f, 1.0f);
		}

		SConfig.Set2DAreaFov(DegreeToRad(85.0f));
	}
	
	~OculusRiftImplemetation() 
	{
		Uninit();
	}
	
	void Uninit()
	{
		// No OVR functions involving memory are allowed after this.
		if(System::IsInitialized())
		{
			RemoveHandlerFromDevices();
			if(nullptr!=pSensor)
			{
				pSensor->Release();
				//delete pSensor; pSensor = nullptr;
			}
			if(nullptr!=pHMD)
			{
				pHMD->Release();
				//delete pHMD; pHMD = nullptr;
			}
			if(nullptr!=pManager)
			{
				pManager->Release();
				//delete pManager; pManager = nullptr;
			}
			System::Destroy();
		}
	}
	
	void OnUpdate()
	{
		// Handle Sensor motion.
		// We extract Yaw, Pitch, Roll instead of directly using the orientation
		// to allow "additional" yaw manipulation with mouse/controller.
		if(pSensor)
		{
			Quatf    hmdOrient = SFusion.GetPredictedOrientation();

			hmdOrient.GetEulerAngles<Axis_Y, Axis_X, Axis_Z>(&m_yaw, &m_pitch, &m_roll);

			/*Player.EyeYaw += (yaw - Player.LastSensorYaw);
			Player.LastSensorYaw = yaw;*/

			// NOTE: We can get a matrix from orientation as follows:
			// Matrix4f hmdMat(hmdOrient);

			// Test logic - assign quaternion result directly to view:
			// Quatf hmdOrient = SFusion.GetOrientation();
			// View = Matrix4f(hmdOrient.Inverted()) * Matrix4f::Translation(-EyePos);
		}
	}

	bool HasSensor() const { return (nullptr!=pSensor); }

	float Yaw() const { return m_yaw; }
	float Pitch() const { return m_pitch; }
	float Roll() const { return m_roll; }

	void GetDistortionValues(float &XCenterOffset, float &Scale, float &K0, float &K1, float &K2, float &K3)
	{
		const DistortionConfig& Distortion = SConfig.GetDistortionConfig();
		//StereoEyeParams eyeParams = SConfig.GetEyeRenderParams();
		XCenterOffset	= Distortion.XCenterOffset;
		Scale			= Distortion.Scale;
		K0				= Distortion.K[0];
		K1				= Distortion.K[1];
		K2				= Distortion.K[2];
		K3				= Distortion.K[3];
	}

	matrix4x4f GetPerspectiveMatrix(const ViewEye eye)
	{
		const float ProjectionCenterOffset = SConfig.GetProjectionCenterOffset();
		const float YFov = SConfig.GetYFOVRadians();
		const float Aspect = SConfig.GetAspect();

		// Projection matrix for the center eye, which the left/right matrices are based on.
		float znear, zfar;
		Pi::renderer->GetNearFarRange(znear, zfar);
		matrix4x4f projCenter = PerspectiveRH(YFov, Aspect, znear, zfar);
		matrix4x4f projLeft  = matrix4x4f::Translation(ProjectionCenterOffset, 0, 0) * projCenter,
				   projRight = matrix4x4f::Translation(-ProjectionCenterOffset, 0, 0) * projCenter;

		matrix4x4f tMat = projCenter;
		switch(eye) {
		case ViewEye_Centre:
		default:
			return projCenter;
		case ViewEye_Left: return projLeft;
		case ViewEye_Right:return projRight;
		}
	}

	float GetYFOVDegrees() { return SConfig.GetYFOVDegrees(); }

	const OVR::HMDInfo& GetHMDInfo() const { return HMDInfo; }

private:
	matrix4x4f PerspectiveRH(const float yfov, const float aspect, const float znear, const float zfar)
	{
		matrix4x4f m(0.0f);
		const float tanHalfFov = tan(yfov * 0.5f);
		m[0] = 1.0f / (aspect * tanHalfFov);
		m[5] = 1.0f / tanHalfFov;
		m[10] = zfar / (znear - zfar);
		m[11] = -1.0f;
		m[14] = (zfar * znear) / (znear - zfar);
		m[15] = 0.0f;     
		return m;
	}

	// *** Oculus HMD Variables
    DeviceManager	*pManager;
    SensorDevice	*pSensor;
    HMDDevice		*pHMD;
    SensorFusion				SFusion;
    OVR::HMDInfo				HMDInfo;

	// Stereo view parameters.
    StereoConfig        SConfig;

	int Width;
    int Height;

	// 
	float	m_yaw;
	float	m_pitch;
	float	m_roll;
};

// static members for the interface class.
std::unique_ptr<OculusRiftImplemetation> OculusRiftInterface::mPimpl;

void OculusRiftInterface::Init()
{
	mPimpl.reset(new OculusRiftImplemetation);
}

void OculusRiftInterface::Uninit()
{
	assert(mPimpl.get());
	mPimpl->Uninit();
}

void OculusRiftInterface::Update()
{
	assert(mPimpl.get());
	mPimpl->OnUpdate();
}

//static
bool OculusRiftInterface::HasHMD()
{
	assert(mPimpl.get());
	return mPimpl->HasSensor();
}

//static 
void OculusRiftInterface::GetYawPitchRoll(float &yaw, float &pitch, float &roll)
{
	assert(mPimpl.get());
	yaw		= mPimpl->Yaw();
	pitch	= mPimpl->Pitch();
	roll	= mPimpl->Roll();
}

//static
void OculusRiftInterface::GetDistortionValues(float &XCenterOffset, float &Scale, float &K0, float &K1, float &K2, float &K3)
{
	assert(mPimpl.get());
	mPimpl->GetDistortionValues(XCenterOffset, Scale, K0, K1, K2, K3);
}

//static 
matrix4x4f OculusRiftInterface::GetPerspectiveMatrix(const ViewEye eye)
{
	assert(mPimpl.get());
	return mPimpl->GetPerspectiveMatrix(eye);
}

//static 
float OculusRiftInterface::GetYFOVDegrees() 
{ 
	return mPimpl->GetYFOVDegrees(); 
}

//static 
OculusRiftInterface::ScreenInfo OculusRiftInterface::GetScreenInfo()
{
	OculusRiftInterface::ScreenInfo si;
	const OVR::HMDInfo& info = mPimpl->GetHMDInfo();
	si.DesktopX = info.DesktopX;
	si.DesktopY = info.DesktopY;
	si.HResolution = info.HResolution;
	si.VResolution = info.VResolution;
	return si;
}
