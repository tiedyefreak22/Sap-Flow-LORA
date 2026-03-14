from Matlab2Python import *
from statistics import mean
from scipy import constants
from scipy.special import erfcinv, erfinv, erf, gammaln, iv
from scipy.spatial.transform import Rotation as R
from scipy.optimize import minimize, fsolve, brentq, newton, bisect, root_scalar
import copy
import pandas as pd
from PD_models import Pd_models_main
import warnings
warnings.filterwarnings("ignore")

# Set the seed for repeatability
np.random.seed(42)

# System parameters
c = constants.c
k_boltzmann = constants.Boltzmann
Rearth = 6371000

class radarDataGenerator:
    """
    radarDataGenerator class.

    Args:
        SensorIndex (int): Unique sensor identifier
        UpdateRate (float): Sensor update rate (Hz)
        DetectionMode (string): Detection mode, specified as "Monostatic", "ESM", or "Bistatic"
        ScanMode (string): Scanning mode of the radar, specified as "Electronic" or "Mechanical"
        TargetReportFormat (string): Format of generated target reports, specified as "Detections", "Tracks", or "Clustered Detections"
        DetectionCoordinates (string): Coordinate system used to report detections, specified as "Scenario", "Body", "Sensor Rectangular", or "Sensor Spherical"
        HasElevation (bool): Enable the radar to scan in elevation and measure target elevation angles
        HasINS (bool): Enable the INS input argument, which passes the current estimate of the sensor platform pose to the sensor
        HasRangeRate (bool): Enable the radar to measure target range rates
        HasRangeAmbiguities (bool): Enable range ambiguities
        HasFalseAlarms (bool): Enable creating false alarm radar measurements
        CenterFrequency (float): Center frequency of the radar band (Hz)
        Bandwidth (float): Radar waveform bandwidth (Hz)
        RangeResolution (float): Range resolution of the radar (m)
        AzimuthResolution (float): Current effective azimuthal resolution of the sensor
        ElevationResolution (float): Current effective elevation resolution of the sensor
        ReferenceRange (float): Reference range for the given probability of detection and the given reference radar cross-section (RCS) (m)
        ReferenceRCS (float): Reference radar cross-section (RCS) for a given probability of detection and reference range (dBsm)
        DetectionProbability (float): Probability of detecting a reference target, specified as a scalar in the range (0, 1]
        FalseAlarmRate (float): False alarm report rate within each radar resolution cell, specified as a positive real scalar in the range [10–7, 10–3] (dimensionless)
        RangeLimits (np.array): 1 x 2 vector in m (e.g. [min, max])
        MaxUnambiguousRange (float): Maximum unambiguous detection range (m)
        FieldOfView (np.array): 2 x 1 vector in degrees (e.g. [[2], [6]])
        MountingAngles (np.array): Mounting rotation angles of the radar, in degrees, specified as 1 x 3 pointing vector e.g. [gamma, beta, alpha]
        ElectronicAzimuthLimits (np.array): 1 x 2 vector in deg (e.g. [min, max])
        ElectronicElevationLimits (np.array): 1 x 2 vector in deg (e.g. [min, max])
        MechanicalAzimuthLimits (np.array): 1 x 2 vector in deg (e.g. [min, max])
        MechanicalElevationLimits (np.array): 1 x 2 vector in deg (e.g. [min, max])
    """
    def __init__(
        self,
        SensorIndex,
        UpdateRate,
        DetectionMode,
        ScanMode,
        TargetReportFormat,
        DetectionCoordinates,
        HasElevation,
        HasINS,
        HasRangeRate,
        HasRangeAmbiguities,
        HasFalseAlarms,
        CenterFrequency,
        Bandwidth,
        RangeResolution,
        AzimuthResolution,
        ElevationResolution,
        ReferenceRange,
        ReferenceRCS,
        DetectionProbability,
        FalseAlarmRate,
        RangeLimits,
        MaxUnambiguousRange,
        FieldOfView,
        MountingAngles,
        ElectronicAzimuthLimits = None,
        ElectronicElevationLimits = None,
        MechanicalAzimuthLimits = None,
        MechanicalElevationLimits = None,
    ):
        self.SensorIndex = SensorIndex
        self.UpdateRate = UpdateRate
        self.DetectionMode = DetectionMode
        self.ScanMode = ScanMode
        self.TargetReportFormat = TargetReportFormat
        self.DetectionCoordinates = DetectionCoordinates
        self.HasElevation = HasElevation
        self.HasINS = HasINS
        self.HasRangeRate = HasRangeRate
        self.HasRangeAmbiguities = HasRangeAmbiguities
        self.HasFalseAlarms = HasFalseAlarms
        self.CenterFrequency = CenterFrequency
        self.Bandwidth = Bandwidth
        self.RangeResolution = RangeResolution
        self.AzimuthResolution = AzimuthResolution
        self.ElevationResolution = ElevationResolution
        self.ReferenceRange = ReferenceRange
        self.ReferenceRCS = ReferenceRCS
        self.DetectionProbability = DetectionProbability
        self.FalseAlarmRate = FalseAlarmRate
        self.RangeLimits = RangeLimits
        self.MaxUnambiguousRange = MaxUnambiguousRange
        self.ElectronicAzimuthLimits = ElectronicAzimuthLimits
        self.ElectronicElevationLimits = ElectronicElevationLimits
        self.MechanicalAzimuthLimits = MechanicalAzimuthLimits
        self.MechanicalElevationLimits = MechanicalElevationLimits
        self.FieldOfView = FieldOfView
        self.MountingAngles = MountingAngles
        self.PointingBasis = rotate_vector(np.array([1, 0, 0]), *self.MountingAngles)
        self.ReferenceSNR = albersheims(self.DetectionProbability, 1, Pfa = self.FalseAlarmRate)
        self.RadarLoopGain = self.ReferenceSNR + 40 * log10(self.ReferenceRange) - self.ReferenceRCS

        # Check if both are missing
        if self.ScanMode.lower() == "electronic" and (self.ElectronicAzimuthLimits is None or self.ElectronicElevationLimits is None):
            raise ValueError("Please provide ELECTRONIC azimuth and elevation limits.")
        # Check if both are present
        elif self.ScanMode.lower() == "mechanical" and (self.MechanicalAzimuthLimits is None or self.MechanicalElevationLimits is None):
            raise ValueError("Please provide MECHANICAL azimuth and elevation limits.")

    def beamPos(self, simTime = None, frame = None):
        """
        Returns the position of the radar beam at either the prescribed time or frame.
    
        Args:
            simTime (float): Time to end simulation.
            frame (float): Update rate of simulation in hertz.

        Returns:
            lookAng (np.array): 1 x 3 pointing vector e.g. [gamma, beta, alpha] in the sensor-local frame.
        """
        # Check if both are missing
        if simTime is None and frame is None:
            raise TypeError("One of 'simTime' or 'frame' is required.")
        # Check if both are present
        elif simTime is not None and frame is not None:
            raise TypeError("Only one of 'simTime' or 'frame' can be provided.")
        
        # Rest of your function logic goes here
        if simTime is not None:
            frame = floor(simTime / self.UpdateRate)
        
        if self.ScanMode.lower() == "electronic":
            numScanPointsAz = floor(np.diff(self.ElectronicAzimuthLimits)[0] / self.FieldOfView[0][0]) + 1
            numScanPointsEl = floor(np.diff(self.ElectronicElevationLimits)[0] / self.FieldOfView[1][0]) + 1
            numScanPoints = numScanPointsAz * numScanPointsEl
        
            scan = floor(frame / numScanPoints)
            scanPointEl = floor((frame % numScanPoints) / numScanPointsAz)
            scanPointAz = (frame % numScanPoints) % numScanPointsAz
            lookAng = [
                0,
                scanPointEl * self.FieldOfView[1][0] + self.ElectronicElevationLimits[0],
                scanPointAz * self.FieldOfView[0][0] + self.ElectronicAzimuthLimits[0],
            ]
        
        elif self.ScanMode.lower() == "mechanical":
            numScanPointsAz = (2 * floor(np.diff(azLimits)[0] / fov[0][0])) + 1
            numScanPointsEl = (2 * floor(np.diff(elLimits)[0] / fov[1][0])) + 1
            numScanPoints = numScanPointsAz * numScanPointsEl
            
            ElIndexList = list(range(ceil(numScanPointsEl / 2) + 1)) + list(range(ceil(numScanPointsEl / 2) - 1, -1, -1))[:-1] # check this -1 indexing
            AzIndexList = list(range(ceil(numScanPointsAz / 2) + 1)) + list(range(ceil(numScanPointsAz / 2), -1, -1))
            
            scan = floor(frame / ((len(ElIndexList) * len(AzIndexList)) / 4))
            scanPointEl = ElIndexList[floor(frame / (len(AzIndexList) / 2)) % len(ElIndexList)]
            scanPointAz = AzIndexList[frame % len(AzIndexList)]
            lookAng = [
                0,
                scanPointEl * self.FieldOfView[1][0] + self.ElectronicElevationLimits[0],
                scanPointAz * self.FieldOfView[0][0] + self.ElectronicAzimuthLimits[0],
            ]
        
        return lookAng

class radarScenario:
    """
    RadarScenario class.

    Args:
        StopTime (float): Time to end simulation.
        UpdateRate (float): Update rate of simulation in hertz.
        SimulationTime (float): Time to set simulation to.
    """
    def __init__(self, StopTime, UpdateRate, SimulationTime = 0):
        self.StopTime = StopTime
        self.UpdateRate = UpdateRate
        self._SimulationTime = SimulationTime
        self.platformsInScene = []
        self.platformIndices = []

    @property
    def SimulationTime(self):
        """
        Get simulation time

        Returns:
            value (float): Current stored simulation time.
        """
        # return copy.deepcopy(self._SimulationTime)
        return self._SimulationTime
        
    @SimulationTime.setter
    def SimulationTime(self, value):
        """
        Set simulation time
        
        Args:
            value (float): Time to set simulation to.
        """
        # if (value > self.StopTime):
        #     raise ValueError("SimulationTime exceeds StopTime")
        self._SimulationTime = value

    def add_to_scene(self, child_instance, platformIndex):
        """
        Method to register a child instance, or replace an existing one.

        Args:
            child_instance (platform): Sensor or target to register within the scenario.
        """
        if child_instance not in self.platformsInScene: # add new
            self.platformsInScene.append(child_instance)
            if platformIndex is None:
                int_set = set(self.platformIndices)
                current_int = 0
                # Iterate and check for the integer's presence in list
                while current_int in int_set:
                    current_int += 1
                self.platformIndices.append(current_int)
                child_instance.PlatformIndex = current_int
            elif platformIndex not in self.platformIndices:
                self.platformIndices.append(platformIndex)
                child_instance.PlatformIndex = current_int
        else: # replace existing
            list_index = self.platformsInScene.index(child_instance)
            self.platformsInScene[list_index] = child_instance
            if platformIndex not in self.platformIndices:
                self.platformIndices[list_index] = platformIndex
                child_instance.PlatformIndex = current_int
            else:
                int_set = set(self.platformIndices)
                current_int = 0
                # Iterate and check for the integer's presence in list
                while current_int in int_set:
                    current_int += 1
                self.platformIndices[list_index] = current_int
                child_instance.PlatformIndex = current_int
    
    def detect(self):
        """
        Gather and return all detections by all sensors

        Returns:
            detections ([detection, detection,...]): Array of gathered detections in the global frame.
        """
        detections = []

        sensors = (platform for platform in self.platformsInScene if hasattr(platform, "Sensors"))
        targets = (platform for platform in self.platformsInScene if not hasattr(platform, "Sensors"))

        numFrames = int(self.StopTime * self.UpdateRate)
        time = np.linspace(0, numFrames - 1, numFrames).T / self.UpdateRate
        
        # iterate through platforms with sensors
        for sensor in sensors:
            obs_pos = np.array(sensor.Position)
            obs_vec = rotate_vector(sensor.Sensors.PointingBasis, *sensor.Sensors.beamPos(frame = frame))
            obs_fov = np.array(sensor.Sensors.FieldOfView)
            lookAng = sensor.Sensors.beamPos(frame = frame)
    
            for target in targets:
                tgtPosTruth = target.Position + np.outer(time, target.Trajectory.Velocity)
                target_pos = np.array(tgtPosTruth[frame])
                wv = sensor.Sensors.CenterFrequency
                Pt_C = (sensor.Sensors.ReferenceRange ** 4 * sensor.Sensors.ReferenceSNR) / (wv ** 2 * 10 ** (sensor.Sensors.ReferenceRCS / 10))
                
                if self.is_target_in_fov(obs_pos, obs_vec, target_pos, obs_fov):
                    sigma = target.Signatures.Pattern
                    R = mag_vec(target_pos - obs_pos)
                    SNR = (Pt_C * wv ** 2 * sigma) / (R ** 4)
                    
                    detections.append(
                        detection(
                            self.SimulationTime, # Time
                            lookAng, # Measurement
                            zeros((np.size(lookAng), np.size(lookAng))), # MeasurementNoise
                            sensor.Sensors.SensorIndex, # SensorIndex
                            [], # ObjectClassID
                            [], # MeasurementParameters
                            objectAttributes(target._PlatformIndex, SNR), # ObjectAttributes
                        )
                    )

        return detections

    def advance(self):
        """
        Advances simulation by one time step (1 / self.UpdateRate)
        """
        if self.UpdateRate == 0:
            pass
            # sensors = (platform for platform in self.platformsInScene if hasattr(platform, "Sensors"))
            # for sensor in sensors:
                # self.UpdateRate == self.platformsInScene[0].Sensors.UpdateRate
        
        self.SimulationTime += (1 / self.UpdateRate)
        
        if self.SimulationTime >= self.StopTime:
            return 0
        else:
            return 1

    def restart(self):
        """
        Resets simulation time to 0.
        """
        self.SimulationTime = 0

    @staticmethod
    def is_target_in_fov(observer_pos, observer_forward, target_pos, fov_angles_degrees):
        """
        Determines if the target is within the specified fov
    
        Args:
            obs_pos (np.array): 1 x 3 position vector e.g. [x, y, z]
            obs_forward (np.array): 1 x 3 pointing vector e.g. [gamma, beta, alpha]
            tgt_pos (np.array): 1 x 3 position vector e.g. [x, y, z]
            fov_angles_degrees: 2 x 1 fov vector i.e. [[az], [el]]
    
        Returns:
            bool: True if the target is within both the horizontal and vertical FOV, False otherwise.
        """
    
        # Calculate and normalize the direction to the target
        direction_to_target = target_pos - observer_pos
        
        # Special case wherein target and observer are colocated
        if mag_vec(direction_to_target) == 0:
            return True
    
        normalized_direction_to_target = norm_vec(direction_to_target)
        
        # Normalize the observer forward vector (if not already normalized)
        normalized_forward = norm_vec(observer_forward)
    
        # Generate the half angles for comparison
        half_h_fov = fov_angles_degrees[0][0] / 2.0
        half_v_fov = fov_angles_degrees[1][0] / 2.0
        
        # Assume a standard coordinate system:
        # X = Forward
        # Y = Left
        # Z = Up
    
        # Project onto the XY plane for horizontal check
        observer_forward_xy = normalized_forward[[0, 1]]
        target_direction_xy = normalized_direction_to_target[[0, 1]]
        xy_dot = np.dot(observer_forward_xy, target_direction_xy)
        xy_mag_prod = mag_vec(observer_forward_xy) * mag_vec(target_direction_xy)
        angle_h_rad = np.arccos(xy_dot / xy_mag_prod)
        angle_h_deg = np.degrees(angle_h_rad)
    
        # Project onto the XZ plane for vertical check
        observer_forward_xz = normalized_forward[[0, 2]]
        target_direction_xz = normalized_direction_to_target[[0, 2]]
        xz_dot = np.dot(observer_forward_xz, target_direction_xz)
        xz_mag_prod = mag_vec(observer_forward_xz) * mag_vec(target_direction_xz)
        angle_v_rad = np.arccos(xz_dot / xz_mag_prod)
        angle_v_deg = np.degrees(angle_v_rad)
    
        # Check if within the half angles
        # The target is in the FOV if the angle is less than half of the total FOV
        # Note: need to add half fov to account for pointing angle being the bottom right corner of the fov window
        is_in_horizontal_fov = angle_h_deg <= half_h_fov
        is_in_vertical_fov = angle_v_deg <= half_v_fov
    
        return is_in_horizontal_fov and is_in_vertical_fov

class platform:
    """
    Platform class to hold both sensors and targets.

    Args:
        scene (radarScenario): Scenario to register the platform within
        Position (np.array): 1 x 3 position vector e.g. [x, y, z]
        Trajectory (kinematicTrajectory): Trajectory class.
    """
    def __init__(self, scene, Position = None, Trajectory = None, Sensors = None):
        self.scene = scene
        self._PlatformIndex = None
        self.scene.add_to_scene(self, None)
        if Position != None:
            self.Position = Position
        elif Trajectory != None:
            self.Position = Trajectory.Position
        else:
            self.Position = [0, 0, 0]
            
        if Trajectory != None:
            self.Trajectory = Trajectory

        if Sensors != None:
            for Sensor in Sensors:
                self.Sensors.append(Sensor)

    @property
    def Sensors(self):
        """
        Get sensors
        
        Returns:
            Sensor (Sensor): Get sensor.
        """
        return self._Sensors

    @Sensors.setter
    def Sensors(self, Sensor):
        """
        Set sensors and assigns platform index
        
        Args:
            Sensor (Sensor): Set sensor.
        """
        self._Sensors = Sensor
        self.PlatformIndex = Sensor.SensorIndex
        self.scene.add_to_scene(self, self._PlatformIndex)

    @property
    def PlatformIndex(self):
        """
        Get PlatformIndex
        
        Returns:
            _PlatformIndex (int): Get PlatformIndex.
        """
        return self._PlatformIndex

    @PlatformIndex.setter
    def PlatformIndex(self, PlatformIndex):
        """
        Set PlatformIndex
        
        Args:
            PlatformIndex (int): Set PlatformIndex.
        """
        self._PlatformIndex = PlatformIndex

class kinematicTrajectory:
    """
    KinematicTrajectory class.

    Args:
        Position (np.array): 1 x 3 position vector e.g. [x, y, z]
        Velocity (np.array): 1 x 3 velocity vector e.g. [vx, vy, vz]
    """
    def __init__(self, Position, Velocity):
        self.Position = Position
        self.Velocity = Velocity

class rcsSignature:
    """
    RcsSignature class.

    Args:
        Pattern (float): Platform RCS profiles in dBsm
    """
    def __init__(self, Pattern):
        self.Pattern = Pattern

class detection:
    """
    Detection class.

    Args:
        Time (float): 
        Measurement (np.array): N-element vector
        MeasurementNoise (np.array): N x N matrix
        SensorIndex (int): Sensor index number
        ObjectClassID (int): 
        MeasurementParameters (structure array):
        ObjectAttributes (structure array):
    """
    def __init__(
        self,
        Time,
        Measurement,
        MeasurementNoise,
        SensorIndex,
        ObjectClassID,
        MeasurementParameters,
        ObjectAttributes,
    ):
        self.Time = Time
        self.Measurement = Measurement
        self.MeasurementNoise = MeasurementNoise
        self.SensorIndex = SensorIndex
        self.ObjectClassID = ObjectClassID
        self.MeasurementParameters = MeasurementParameters
        self.ObjectAttributes = ObjectAttributes

class objectAttributes:
    """
    ObjectAttributes class.
    
    Args:
        TargetIndex (int): Index of target
        SNR (float): Signal-to-noise ratio
    """
    def __init__(
        self,
        TargetIndex,
        SNR
    ):
        self.TargetIndex = TargetIndex
        self.SNR = SNR

def albersheims(Pd = None, N = None, SNR = None, Pfa = None):
    """
    Calculates the missing value of the provided arguments using Albersheim's equation.

    Args:
        Pd (float): Probability of detection
        N (int): Number of non-coherent pulses
        SNR (float): Signal-to-noise ratio
        Pfa (float): Probability of false alarm

    Returns:
        float: Value not provided above.
    """
    if sum(arg is None for arg in (Pd, N, SNR, Pfa)) == 1:
        if Pfa is None:
            B = log(Pd / (1 - Pd))
            Z = (SNR + 5 * log10(N)) / (6.2 + 4.54 / sqrt(N + 0.44))
            A = (10 ** Z - 1.7 * B) / (1 + 0.12 * B)
            Pfa = 0.62 * exp(-A)
            return Pfa
        elif Pd is None:
            A = log(0.62 / Pfa)
            Z = (SNR + 5 * log10(N)) / (6.2 + 4.54 / sqrt(N + 0.44))
            B = (10 ** Z - A) / (1.7 + 0.12 * A)
            Pd = 1 / (1 + exp(-B))
            return Pd
        elif N is None:
            print("Note: 2 <= N <= 100 in this approximation")
            A = log(0.62 / Pfa)
            B = log(Pd / (1 - Pd))
            Z = log10(A + 0.12 * A * B + 1.7 * B)
            roots = np.roots([5, SNR - 5.5606 * Z, -2.4194 * Z, 0.4125 * Z])
            N = round(max(10 ** roots))
            return N            
        elif SNR is None:
            A = log(0.62 / Pfa)
            B = log(Pd / (1 - Pd))
            SNR = -5 * log10(N) + (6.2 + (4.54 / sqrt(N + 0.44))) * log10(A + 0.12 * A * B + 1.7 * B)
            return SNR
    else:
        raise ValueError("Only one argument allowed to remain 'None'")

def detectability(PDs, PFAs, N = 1, SW = "Swerling0"):
    '''
    Returns the detectability factor of a single radar pulse.
    It leverages the Pd_models library swerling functions (which return Pd),
    but then back solves for the SNR (in dB) that corresponds to that Pd value.

    Args:
        PD (float): 
        PFA (float): 
        N (int): 
        SW (string): 

    Returns:
        np.array: Detectability factor, returned as a J-by-K matrix in dB
    '''
    output_array = []

    if isscalar(PDs):
        PDs = [PDs]
    if isscalar(PFAs):
        PFAs = [PFAs]
    
    for PD in PDs:
        PFA_output_array = []
        for PFA in PFAs:
            Yb = Pd_models_main._calculate_Yb(np.array([N]), np.array([PFA]))
        
            def objective_function(snr_pr):
                match int(SW.lower().replace("swerling", "").strip()):
                    case 1:
                        calc_PD = Pd_models_main._swerling1_calculation(np.array([snr_pr]), np.array([N]), Yb)
                    case 2:
                        calc_PD = Pd_models_main._swerling2_calculation(np.array([snr_pr]), np.array([N]), Yb)
                    case 3:
                        calc_PD = Pd_models_main._swerling3_calculation(np.array([snr_pr]), np.array([N]), Yb)
                    case 4:
                        calc_PD = Pd_models_main._swerling4_calculation(np.array([snr_pr]), np.array([N]), Yb)
                    case 5 | 0:
                        calc_PD = Pd_models_main._swerling0_calculation(np.array([snr_pr]), np.array([N]), Yb)
                    case _:
                        raise ValueError(
                            f"That target type ({target_type}) is not recognized. Please correct it. The form is 'Swerling X' or 'SwerlingX', and capitalization does not matter."
                        )
                return calc_PD - PD

            snrmax = 1e10
            snrmin = 1e-10
            # snr_pr = np.array([1])
            res = brentq(objective_function, snrmin, snrmax)
            output_val = pow2db(res)
            PFA_output_array.append(output_val)
        output_array = np.concatenate((output_array, PFA_output_array), axis=0)
    
    return unlevel(output_array)

def radareqsnr(
    wv,
    tgtrng,
    Pt,
    tau,
    RCS = 1,
    Ts = 290,
    Gain = [20, 20],
    Loss = 0,
    AtmosphericLoss = 0,
    PropagationFactor = 0,
    CustomFactor = 0,
):
    '''
    SNR estimate from radar equation
    
    Args:
        wv (float): Radar wavelength in meters
        tgtrng — two-element row vector of positive values, Range from the receiver to the target. If the radar is monostatic, the transmitter and receiver ranges are identical.
        Pt — Peak transmit power in watts
        tau — Single pulse duration, positive scalar
        RCS — Target's nonfluctuating radar cross section in square meters
        Ts — System noise temperature 290 (default) | positive scalar
        Gain — real-valued 1-by-2 row vector Transmit and Receiver antenna gain. If the radar is monostatic, the transmit and receive antenna gains are identical.
        Loss — General loss factor in decibels that accounts for both system and propagation loss
        AtmosphericLoss — Atmospheric absorption loss 0 (default) | scalar | two-element row vector of real values | length-J column vector of real values | J-by-2 matrix of real values
        PropagationFactor — Propagation factor 0 (default) | scalar | two-element row vector of real values | length-J column vector of real values | J-by-2 matrix of real values
        CustomFactor — Custom factor 0 (default) | scalar | length-J column vector of real values

    Returns:
        SNR (float): SNR in decibels
    '''
    # All of this garbage ensures that if the user only supplied a single value,
    # the necessary data array sizes are filled
    if is_single_row(tgtrng):
        tgtrng = singcolvec(tgtrng)
        
    if is_single_row(PropagationFactor):
        PropagationFactor = singcolvec(PropagationFactor)

    if is_single_row(AtmosphericLoss):
        AtmosphericLoss = singcolvec(AtmosphericLoss)
    
    if (np.size(tgtrng) == 1) or is_single_column(tgtrng):
        tgtrng = np.column_stack((tgtrng, tgtrng))

    if (np.size(PropagationFactor) == 1) or is_single_column(PropagationFactor):
        PropagationFactor = np.column_stack((PropagationFactor, PropagationFactor))

    if (np.size(AtmosphericLoss) == 1) or is_single_column(AtmosphericLoss):
        AtmosphericLoss = np.column_stack((AtmosphericLoss, AtmosphericLoss))

    if (np.size(Gain) == 1) or is_single_column(Gain):
        Gain = unlevel(np.column_stack((Gain, Gain)))
    
    SNR = [(Pt * tau * RCS * wv ** 2) / ((4 * pi) ** 3 * k_boltzmann * Ts * (float(i[0]) * float(i[1])) ** 2) for i in tgtrng]

    if np.shape(SNR)[0] > 1:
        if np.shape(PropagationFactor)[0] > 1 and np.shape(AtmosphericLoss)[0] > 1:
            SNRdb = [pow2db(k) + Gain[0] + Gain[1] + i[0] + i[1] + CustomFactor - Loss - j[0] - j[1] for i, j, k in zip(PropagationFactor, AtmosphericLoss, SNR)]
        elif np.shape(PropagationFactor)[0] > 1 and np.shape(AtmosphericLoss)[0] == 1:
            SNRdb = [pow2db(k) + Gain[0] + Gain[1] + i[0] + i[1] + CustomFactor - Loss - AtmosphericLoss[0][0] - AtmosphericLoss[0][1] for i, k in zip(PropagationFactor, SNR)]
        elif np.shape(PropagationFactor)[0] == 1 and np.shape(AtmosphericLoss)[0] > 1:
            SNRdb = [pow2db(k) + Gain[0] + Gain[1] + PropagationFactor[0][0] + PropagationFactor[0][1] + CustomFactor - Loss - j[0] - j[1] for j, k in zip(AtmosphericLoss, SNR)]
        else:
            SNRdb = [pow2db(k) + Gain[0] + Gain[1] + PropagationFactor[0][0] + PropagationFactor[0][1] + CustomFactor - Loss - AtmosphericLoss[0][0] - AtmosphericLoss[0][1] for k in SNR]
    else:
        if np.shape(PropagationFactor)[0] > 1 and np.shape(AtmosphericLoss)[0] > 1:
            SNRdb = [pow2db(SNR) + Gain[0] + Gain[1] + i[0] + i[1] + CustomFactor - Loss - j[0] - j[1] for i, j in zip(PropagationFactor, AtmosphericLoss)]
        elif np.shape(PropagationFactor)[0] > 1 and np.shape(AtmosphericLoss)[0] == 1:
            SNRdb = [pow2db(SNR) + Gain[0] + Gain[1] + i[0] + i[1] + CustomFactor - Loss - AtmosphericLoss[0][0] - AtmosphericLoss[0][1] for i in PropagationFactor]
        elif np.shape(PropagationFactor)[0] == 1 and np.shape(AtmosphericLoss)[0] > 1:
            SNRdb = [pow2db(SNR) + Gain[0] + Gain[1] + PropagationFactor[0][0] + PropagationFactor[0][1] + CustomFactor - Loss - j[0] - j[1] for j in AtmosphericLoss]
        else:
            SNRdb = pow2db(SNR) + Gain[0] + Gain[1] + PropagationFactor[0][0] + PropagationFactor[0][1] + CustomFactor - Loss - AtmosphericLoss[0][0] - AtmosphericLoss[0][1]
        
    return SNRdb

def beamloss(is2d = False):
    '''
    Beam shape loss for Gaussian antenna pattern

    Args:
        is2d (bool): 

    Returns:
        float: beam shape loss for a radar that scans over one or two angular dimension(s) (dB)
    '''
    if not is2d:
        return pow2db(sqrt((8 * log(2)) / pi))
    else:
        return pow2db((8 * log(2)) / pi)

def arrayscanloss(PD, PFA, N, THETAM = [-60, 60], SW = "Swerling0", COSINEPOWER = 2.5):
    '''
    Loss due to electronic scanning off broadside

    Args:
        PD (np.array): Desired probability of detection, specified as a scalar or J-length vector between 0.1 and 0.999999
        PFA (np.array): Probability of false alarm, specified as a scalar or K-length vector between 1e-15 and 1e-3
        N: (int): Number of received pulses
        THETAM (np.array): Scan sector limits, specified as a scalar or two-element vector
        SW (string): Scan sector limits, specified as the Swerling case for the chi-square distributed target
        COSINEPOWER (float): Exponent of the cosine modeling the gain loss of an array scanned off broadside

    Returns:
        LSS (np.array): Two-way statistical scan sector loss, returned as a J-by-K matrix, where J and K are the dimensions of the PD and PFA arguments. Units are in decibels (dB).
    '''    
    # % The statistical scan sector loss is computed by averaging Pd over the SNR
    # % weighted by the reduction in the array gain due to scanning off-broadside
    
    # % Size of the scan sector in degrees
    dThetam = abs(THETAM[1] - THETAM[0])
    
    # % Sample the scan sector at one degree intervals or take total 10 samples
    # % if the size of the sector is less than 10 degrees
    M = max(dThetam, 10)
    theta = np.linspace(THETAM[0], THETAM[1], M)
    
    # % Reduction in the gain of an array scanned off-broadside by the angle theta
    H = cos(np.radians(theta)) ** COSINEPOWER
    
    # % Statistical scan sector loss
    Lss = statisticalloss(PD, PFA, N, SW, theta, H)

    # % To support arbitrary input dimensions in codegen
    # Lss_col = np.ravel(Lss)
    # lss_idxs = np.isfinite(Lss)
    # Lss[lss_idxs] = pow2db(Lss_col[lss_idxs])

    return Lss

def statisticalloss(Pd, Pfa, N, SW, X, Y):
    '''
    Computes the statistical loss L (in dB) defined as an increase in the input SNR required to maintain
    the desired probability of detection PD, as PD varies as a function of a variable X representing the
    target position in one of the four dimensions of the radar space.

    Args:
        PD (np.array): Desired probability of detection, specified as a scalar or J-length vector between 0.1 and 0.999999
        PFA (np.array): Probability of false alarm, specified as a scalar or K-length vector between 1e-15 and 1e-3
        N: (int): Number of received pulses (must be positive)
        SW (string): Scan sector limits, specified as the Swerling case for the chi-square distributed target
        X (np.array): length-M vector such that X is a grid of values sampled from one of the four dimensions of the radar space
        Y (np.array): length-M vectors such that Y is a corresponding power ratio of the output signal to a reference value that would apply if the sensitivity to X was absent.

    Returns:
        L (np.array): The output is a J-by-K matrix with rows corresponding to PD and columns corresponding to PFA.
    '''  
    if np.size(Pd) > 1:
        Pd = np.array(Pd)[:]
    elif np.size(Pd) == 1:
        Pd = np.array(Pd)
    else:
        raise ValueError("Pd does not contain values")

    if np.size(Pfa) > 1:
        Pfa = np.array(Pfa)[:].T
    elif np.size(Pfa) == 1:
        Pfa = np.array(Pfa)
    else:
        raise ValueError("Pfa does not contain values")
    
    # % The number of available independent target samples for Swerling case SW
    # % and the number of received pulses N
    Ne = swerlingdof(N, SW)
    
    # % Detectability when there is no sensitivity to X
    D = bartonSNRLinear(Pd, Pfa, N, Ne)
    
    # % The lower and the upper bound on the SNR when the statistical loss is
    # % present
    snrLowerBound = D / max(Y)
    snrUpperBound = D / min(Y)
    
    numelPd = np.size(Pd)
    numelPfa = np.size(Pfa)
    Ds = zeros((numelPd, numelPfa))
    dX = abs(X[-1] - X[0])

    # % Solve for the required SNR for each combination of PD and PFA
    if (numelPd > 1) and (numelPfa > 1):
        for i in range(0, numelPd):
            for j in range(0, numelPfa):
                # % Find the SNR such that the average PD computed for SNR*Y(X)
                # % equals to the desired PD
                # % For numerical stability it is important that Y(X)~=0
                def func(snr):
                    return np.trapz(bartonPdLinear(snr * Y, Pfa[j], N, Ne), X) / dX - Pd[i]
        
                if np.isfinite(snrLowerBound[i, j]) and np.isfinite(snrUpperBound[i, j]):
                    Ds[i, j] = brentq(func, snrLowerBound[i, j], snrUpperBound[i, j])

    # Something is up with this block
    elif (numelPd == 1) and (numelPfa > 1):
        Ds = []
        for i in range(0, numelPfa):
            # % Find the SNR such that the average PD computed for SNR*Y(X)
            # % equals to the desired PD
            # % For numerical stability it is important that Y(X)~=0
            def func(snr):
                return np.trapz(bartonPdLinear(snr * Y, Pfa[i], N, Ne), X) / dX - Pd
    
            if np.isfinite(snrLowerBound[i]) and np.isfinite(snrUpperBound[i]):
                Ds.append(brentq(func, snrLowerBound[i], snrUpperBound[i]))
    elif (numelPd > 1) and (numelPfa == 1):
        Ds = []
        for i in range(0, numelPd):
            # % Find the SNR such that the average PD computed for SNR*Y(X)
            # % equals to the desired PD
            # % For numerical stability it is important that Y(X)~=0
            def func(snr):
                return np.trapz(bartonPdLinear(snr * Y, Pfa, N, Ne), X) / dX - Pd[i]

            if np.isfinite(snrLowerBound[i]) and np.isfinite(snrUpperBound[i]):
                Ds.append(brentq(func, snrLowerBound[i], snrUpperBound[i]))
    else:
        # % Find the SNR such that the average PD computed for SNR*Y(X)
        # % equals to the desired PD
        # % For numerical stability it is important that Y(X)~=0
        def func(snr):
            return np.trapz(bartonPdLinear(snr * Y, Pfa, N, Ne), X) / dX - Pd

        if np.isfinite(snrLowerBound) and np.isfinite(snrUpperBound):
            Ds = brentq(func, snrLowerBound, snrUpperBound)
    # % Loss
    L = Ds / D
    return L

def swerlingdof(N, SwerlingCase):
    '''
    Returns the number of degrees of freedom in the chi-squared distributed radar echo received from a Swerling target.

    Args:
        is2d (bool): 
        N (np.array): length-J vector indicating the number of received pulses
        SwCase (string): can take one of the following values 'Swerling0' | 'Swerling1' | 'Swerling2' | 'Swerling3' | 'Swerling4' | 'Swerling5'

    Returns:
        Ne (np.array): 
    '''
    match int(SwerlingCase.lower().replace("swerling", "").strip()):
        case 1:
            Ne = ones(np.shape(N))
        case 2:
            Ne = N
        case 3:
            Ne = 2 * ones(np.shape(N))
        case 4:
            Ne = 2 * N
        case 5 | 0:
            Ne = np.full(np.shape(N), inf)
        case _:
            raise ValueError(
                f"That target type ({target_type}) is not recognized. Please correct it. The form is 'Swerling X' or 'SwerlingX', and capitalization does not matter."
            )
    return Ne

def bartonSNRLinear(Pd, Pfa, N, Ne):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %bartonSNRLinear Barton's Universal Equation for SNR
    %   SNR = bartonSNRLinear(PD,PFA,N,NE) returns the signal-to-noise ratio,
    %   SNR (in linear units), for a given probability of detection, PD, the
    %   probability of false alarm, PFA, the number of received pulses, N, and
    %   a number of independent samples, NE, received from a chi-squared
    %   distributed target with 2*NE degrees of freedom assuming a square law
    %   detector. The input PD is a scalar or a length-J vector, PFA is a
    %   scalar or a length-K vector, and N and NE are scalars. The output SNR
    %   is a J-by-K matrix with rows corresponding to PD and columns
    %   corresponding to PFA.
    %
    %   % Example: Compute the required SNR as a function of the detection
    %   %   probability given the probability of false alarm of 1e-6 and 24
    %   %   pulses received from a Swerling 2 case target.
    %
    %   Pfa = 1e-6                 % Probability of false alarm
    %   Pd = 0.01 : 0.01 : 0.99    % Probability of detection
    %   N = 24                     % Number of received pulses
    %
    %   SNR = bartonSNRLinear( Pd, Pfa, N, 2 )
    
    %#codegen
    '''
    if np.size(Pd) > 1:
        Pd = np.array(Pd)[:]
    elif np.size(Pd) == 1:
        Pd = np.array(Pd)
    else:
        raise ValueError("Pd does not contain values")

    if np.size(Pfa) > 1:
        Pfa = np.array(Pfa)[:].T
    elif np.size(Pfa) == 1:
        Pfa = np.array(Pfa)
    else:
        raise ValueError("Pfa does not contain values")

    if np.isfinite(Ne):
        Tpfa = real(gammaincinv(1 - Pfa, N))
        Tpd = real(gammaincinv(1 - Pd, Ne))
        SNR = ((1 / Tpd) * (Tpfa - (N - Ne)) - 1) * (Ne / N)
    else:
        # % Use Shnidman's approximation for a non-fluctuating target
        SNR = shnidmanNonfluctuating(Pd, Pfa, N)

    SNR = np.maximum(SNR, 0)
    return SNR

def shnidmanNonfluctuating(Pd, Pfa, N):
    '''
    %This function is for internal use only. It may be removed in the future.
    % shnidmanNonfluctuating is called from shnidman in phased array toolbox
    % and from bartonsnr in radar toolbox.
    % When this function is called from bartonsnr, Pd is a column vector and
    % Pfa is a row vector. The output is expected to be a numel(Pd) by
    % numel(Pfa) matrix.
    
    %shnidmanNonfluctuating shnidman's approximation for Swerling 0/5 case
    %   SNR = shnidmanNonfluctuating(PD,PFA,N) returns the signal-to-noise
    %   ratio, SNR (in linear units), for a given probability of detection, PD,
    %   the probability of false alarm, PFA, and the number of received pulses,
    %   N, assuming a nonfluctuating target. The input PD is a scalar or a
    %   length-J vector, PFA is a scalar or a length-K vector, and N is a
    %   scalar. The output SNR is a J-by-K matrix with rows corresponding to PD
    %   and columns corresponding to PFA.
    %
    %   % Example: Compute the required SNR as a function of the detection
    %   %   probability given the probability of false alarm of 1e-6 and 24 
    %   %   received pulses.
    %
    %   Pfa = 1e-6                 % Probability of false alarm
    %   Pd = 0.01 : 0.01 : 0.99    % Probability of detection
    %   N = 24                     % Number of received pulses
    %
    %   SNR = shnidmanNonfluctuating( Pd, Pfa, N )
    
    %   Copyright 2020 The MathWorks, Inc.
    
    %   Reference
    %   [1] Mark Richards, Fundamentals of Radar Signal Processing, Page 337
    
    %#codegen
    '''
    alpha = 0
    if N > 40:
        alpha = 0.25
    elif N < 0.5:
        alpha = 0.25 - N / 2

    etaPfa = sqrt(-0.8 * log(4 * Pfa * (1 - Pfa)))
    etaPd = sign(Pd - 0.5) * sqrt(-0.8 * log(4 * Pd * (1 - Pd)))

    # % When called from bartonsnr, Pd is a column vector and Pfa is a row vector
    if isscalar(Pd) and isscalar(Pfa):
        eta = etaPd + etaPfa
    else:
        eta = np.tile(etaPd, np.shape(Pfa)) + np.tile(etaPfa, np.shape(Pd))

    SNR = (eta / N) * (eta + 2 * sqrt(N / 2 + (alpha - 0.25)))
    SNR = np.maximum(SNR, 0)

    return SNR

def bartonPdLinear(SNR, Pfa, N, Ne):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %bartonPdLinear Barton's Universal Equation for probability of detection
    %   PD = bartonPdLinear(SNR,PFA,N,NE) returns the probability of detection,
    %   PD, given the signal-to-noise ratio, SNR (in linear units), the
    %   probability of false alarm, PFA, the number of received pulses, N, and
    %   a number of independent samples, NE, received from a chi-squared
    %   distributed target with 2*NE degrees of freedom assuming a square law
    %   detector. The input SNR is a scalar or a length-J vector, PFA is a
    %   scalar or a length-K vector, and N and NE are scalars. The output PD is
    %   a J-by-K matrix with rows corresponding to SNR and columns
    %   corresponding to PFA.
    %
    %   % Example: Compute the detection probability as a function of the
    %   %   signal-to-noise ratio given the false alarm probability of 1e-6
    %   %   and 24 pulses received from a Swerling 2 case target.
    %
    %   Pfa = 1e-6           % Probability of false alarm
    %   N = 24               % Number of received pulses
    %   SNR = -5 : 0.5 : 10  % dB
    %
    %   Pd = bartonPdLinear( SNR, Pfa, N, 2 )
    
    %#codegen
    '''
    if np.size(SNR) > 1:
        SNR = np.array(SNR)[:]
    elif np.size(SNR) == 1:
        SNR = np.array(SNR)
    else:
        raise ValueError("SNR does not contain values")

    if np.size(Pfa) > 1:
        Pfa = np.array(Pfa)[:].T
    elif np.size(Pfa) == 1:
        Pfa = np.array(Pfa)
    else:
        raise ValueError("Pfa does not contain values")
    
    if np.isfinite(Ne):
        T = real(gammaincinv(N, 1 - Pfa))
        Pd = 1 - real(gammainc((1 / (SNR * N / Ne + 1)) * (T - (N - Ne)), Ne))
    else:
        # % Use Shnidman's approximation for a non-fluctuating target
        alpha = 0
        if N > 40:
            alpha = 0.25
        elif N < 0.5:
            N = 0.5
    
        etaPfa = sign(0.5 - Pfa) * sqrt(-0.8 * log(4 * Pfa * (1 - Pfa)))
        etaPd = sqrt(N / 2 + alpha - 0.25) - sqrt(SNR * N + N / 2 + alpha - 0.25)
        eta = etaPfa + etaPd
    
        Pd = 0.5 * (1 - sign(eta) * sqrt(1 - exp(-5 / 4 * (eta ** 2))))
    
    Pd[Pd > 1 ] = 1
    Pd[Pd < 0] = 0
    return Pd

def radarpropfactor(
    R,
    freq,
    anht,
    tgtht,
    pol = "H",
    epsc = None,
    htsd = 0.01,
    beta0 = None,
    vegTypeP = "None",
    hpbw = 10,
    anpat = None,
    patel = None,
    tilt = 0,
    Re = None,
    n0 = 1.000318
):
    '''
    pol = Polarization
    epsc = SurfaceRelativePermittivity
    htsd = SurfaceHeightStandardDeviation
    beta0 = SurfaceSlope
    hpbw = ElevationBeamwidth
    anpat = AntennaPattern
    patel = PatternAngles
    tilt = TiltAngle
    vegTypeP = VegetationType
    Re = EffectiveEarthRadius
    n0 = RefractiveIndex
    
    radarpropfactor One-way radar propagation factor 
       F = radarpropfactor(R,FREQ,ANHT) calculates the one-way propagation
       factor F in dB assuming a surface target and a sea state of 0. The
       propagation takes into account the interference between the direct and
       ground-reflected rays.
    
       The calculation is performed for an M-length vector of free space
       ranges R in meters at the positive, scalar operating frequency FREQ in
       Hz for a scalar antenna height ANHT in meters. The propagation factor F
       is an M-length column vector in dB. The heights are assumed to be
       referenced to the surface.
    
       The calculation estimates the complex relative permittivity (dielectric
       constant) of the reflecting surface using a sea water model that is
       valid from 100 MHz to 10 GHz.
     
       The target height is assumed to be the height of significant clutter
       sources above the average surface height. Specifically, the target
       height is calculated as 3 times the standard deviation of the surface
       height.
    
       Atmospheric refraction is taken into account through the use of an
       effective Earth radius. Scattering and ducting are assumed to be
       negligible.
     
       Assuming the paths are the same, the two-way propagation factor can be
       calculated by multiplying by 2 the one-way propagation factor output by
       this function.
     
       F = radarpropfactor(R,FREQ,ANHT,TGTHT) calculates the target
       propagation factor F in dB assuming a scalar or M-length target height
       of TGTHT in meters.
    
       F = radarpropfactor(...,'Polarization',POL) specifies the polarization
       of the transmitted wave. The parameter, POL, can be one of 'H' | 'V'.
       'H' indicates horizontal polarization and 'V' indicates vertical
       polarization. POL defaults to 'H'.
     
       F = radarpropfactor(...,'SurfaceRelativePermittivity',EPSC) specifies
       the scalar, complex relative permittivity (dielectric constant) of the
       reflecting surface. The default value of EPSC depends on the value of
       FREQ. The function uses a sea water model that is valid up to 10 GHz.
     
       F = radarpropfactor(...,'SurfaceHeightStandardDeviation',HGTSD)
       specifies the scalar standard deviation of the surface height in
       meters. The default value of HGTSD is 0.01 meters, indicating a sea
       state of 0.
    
       F = radarpropfactor(...,'SurfaceSlope',BETA0) specifies the nonnegative
       scalar surface slope in degrees. This value is expected to be 1.4 times
       the RMS surface slope. Given the condition that
            2*GRAZ/BETA0 < 1, 
       where GRAZ is the grazing angle of the geometry specified in degrees,
       the effective surface height standard deviation in meters is calculated
       as
            Effective HGTSD = HGTSD*(2*GRAZ/BETA0)^(0.2), 
       which better accounts for shadowing. Otherwise, the effective height
       standard deviation is equal to HGTSD. BETA0 defaults to the surface
       slope value output by the searoughness function for a sea state of 0.
       
       F = radarpropfactor(...,'VegetationType',VEGTYPE) specifies the
       vegetation type of the surface as one of either 'Trees' | 'Weeds' |
       'Brush' | 'Grass' | 'None'. In the case of VEGTYPE set to 'Trees' |
       'Weeds' | 'Brush', there is an assumption of dense vegetation. The case
       of VEGTYPE set to 'Grass' assumes thin grass. Use this argument when
       using the function on surfaces different from the sea. Defaults to
       'None'.
     
       F = radarpropfactor(...,'ElevationBeamwidth',ELBW) specifies the
       positive scalar half-power elevation beamwidth in degrees. The
       elevation beamwidth is used in the calculation of a sinc antenna
       pattern. The default antenna pattern is symmetrical with respect to the
       beam maximum and is of the form sin(u)/u. The parameter u is given by
            u = k*sin(theta), 
       where theta is the elevation angle in radians and k is given by
            k = 1.39157/sin(ELBW/2). 
       ELBW defaults to 10 degrees.
    
       F = radarpropfactor(...,'AntennaPattern',PAT,'PatternAngles', PATEL)
       specifies the antenna elevation pattern's normalized voltage response
       in linear units (Volts) and corresponding elevation angles in degrees.
       This is an alternative to specifying the elevation beamwidth. Both PAT
       and PATEL must be vectors of the same size. PATEL must be between -90
       and 90 degrees. In general, to properly compute the coverage, the
       pattern should be specified from -90 to 90 degrees. If both an antenna
       pattern and an elevation beamwidth are provided, the function uses the
       antenna pattern and ignores the elevation beamwidth value. Defaults to
       a sinc antenna pattern.
     
       F = radarpropfactor(...,'TiltAngle',TILTANG) specifies the scalar tilt
       angle in degrees of the antenna with respect to the surface. TILTANG
       defaults to 0.
    
       F = radarpropfactor(...,'EffectiveEarthRadius',RE) specifies the
       effective Earth radius as a positive scalar in meters. The effective
       Earth radius is an approximation used for modeling refraction effects
       in the troposphere. The default calculates the effective Earth radius
       using a refractivity gradient of -39e-9, which results in approximately
       4/3 of the real Earth radius.
    
       F = radarpropfactor(...,'RefractiveIndex',N0) specifies the refractive
       index at the surface as a nonnegative scalar. Defaults to approximately
       1.000318, which is the output of the refractiveidx function at an
       altitude of 0 meters.
    
       radarpropfactor(...) plots the one-way propagation factor (dB) versus
       range in km. Default range units are km.
    
        Examples:
    
        Example 1:
          Plot the propagation factor for a surface target as seen by a 10 
          GHz X-band radar assuming an antenna height of 20 m.
       R     = (0.1:0.01:100).*1e3 % Range (m)
       freq  = 10e9                % Frequency (Hz) 
       anht  = 20                  % Radar height (m)
       radarpropfactor(R,freq,anht)
    
        Example 2:
          Plot the propagation factor for a 3 GHz S-band radar assuming an
          antenna height of 10 m and a target height of 1 km. Assume that the
          surface has a height standard deviation of 1 m, and the surface 
          slope is 0.05 deg. 
       R     = (30:0.5:180)*1e3 % Range (m)
       freq  = 3e9              % Frequency (Hz) 
       anht  = 10               % Radar height (m)
       tgtht = 1e3              % Target height (m) 
       hgtsd = 1                % Height standard deviation (m) 
       beta0 = 0.05             % Surface slope (deg)
       radarpropfactor(R,freq,anht,tgtht,...
           'SurfaceHeightStandardDeviation',hgtsd,...
           'SurfaceSlope',beta0)
     
       See also searoughness, landroughness, el2height, radareqsnr, radarvcd,
       blakechart, earthSurfacePermittivity, buildingMaterialPermittivity.
    
       Copyright 2020-2023 The MathWorks, Inc.
    
       References
         [1] Blake, L.V. "Machine Plotting of Radar Vertical-Plane Coverage
             Diagrams." Naval Research Laboratory, 1970 (NRL Report 7098).
         [2] Barton, David K. Radar Equations for Modern Radar. Norwood, MA:
             Artech House, 2013.
    '''
    if epsc is None:
        epsc = seaComplexPermittivity(freq)
    if beta0 is None:
        _, beta0 = searoughness(0)
    if (anpat is None) and (patel is None):
        anpat, patel = sincpattern(hpbw)
    elif (anpat is None):
        anpat, _ = sincpattern(hpbw)
    else:
        _, patel = sincpattern(hpbw)
    if Re is None:
        Re, _ = effearthradius()
        
    if isscalar(tgtht):
        tgtht = [tgtht]

    rfs = np.sort(np.ravel(R))
    numTgt = numel(tgtht)
    numR = numel(rfs)
    FdiffSamp = ones((numR, 1))
    Fsamp = ones((numR, 1))
    FdB = zeros(numR)
    for it in range(0, numTgt):
        tgtht = tgtht[it]
        thisRfs, idx = radarpropfactor_getRange(rfs, it, numTgt)
    
        # Calculate range for start of diffraction region
        Rhoriz = sqrt(2 * Re) * (sqrt(anht) + sqrt(tgtht)) # Start of diffraction region, Eqn 41, ref 1
    
        # Calculate composite propagation factor
        match pol:
            # H-pol propagation factor composed of:
            #    1. Interference region
            #    2. Intermediate region (if thisRfs(end) > Rdelta)
            #    3. Diffraction region (if thisRfs(end) >= Rhoriz)
            case 'H':
                # Over-sample range
                matlabmax(np.hstack((thisRfs, Rhoriz)), [], 1)
                rmax = matlabmax(singcolvec(np.hstack((thisRfs, Rhoriz))), [], 1)[0]
                rmin = np.absolute(tgtht - anht)
                rngSamp = radarpropfactor_getRangeSamples(thisRfs, rmax, rmin)
    
                # Calculate direct ray elevation angles
                thetad, rngiSamp, maxHt, minHt, flipGeom = directRayElAng(rngSamp, anht, tgtht, Rhoriz, Re)
    
                # Shift antenna pattern
                anpat, patel, tiltOut = antennaPatternShift(anpat, patel, tilt, flipGeom)
    
                # Recalculate interference region using range samples
                FinterfSamp, pdiSamp, wv = propfactorinterf(freq, minHt, thetad, epsc, pol, anpat, patel, tiltOut, htsd, beta0, vegTypeP, Re)
    
                # Calculate Rdelta
                if tgtht <= 3 * htsd:
                    # Surface target calculation (approximation)
                    Rdelta = Re * wv / (12 * tgtht) * (sqrt(1 + (2 * anht) / Re * ((12 * tgtht / wv) ** 2)) - 1) # Eqn 9.12, ref 2
                else:
                    # Interpolate to find point = wv/6
                    x = wv / 6
                    test_var, idxNearest = matlabmin(np.absolute(pdiSamp - x), [], 2)
                    # idxInterp = range(max(idxNearest - 1, 1), min(idxNearest + 1, numel(pdiSamp)))
                    idxInterp = np.arange(np.maximum(idxNearest - 1, 0), np.minimum(idxNearest + 1, len(pdiSamp) - 1) + 1)
                    ymax = rngiSamp[idxInterp[-1]]
                    ymin = rngiSamp[idxInterp[0]]
                    xmax = pdiSamp[idxInterp[-1]]
                    xmin = pdiSamp[idxInterp[0]]
                    Rdelta = ymin + (ymax - ymin) * (x - xmin) / (xmax - xmin)

                # Calculation of diffraction and intermediate regions
                if thisRfs[-1] > Rdelta:
                    # Calculate diffraction region
                    FdiffSamp = propfactordiff(rngSamp, freq, minHt, maxHt, Re, n0)
    
                    # Calculate interference + intermediate + diffraction
                    Fsamp = propfactorall(rngSamp, FinterfSamp, FdiffSamp, Rdelta, Rhoriz)
    
                    # Interpolate values to user requested ranges
                    F = interp1(rngSamp, Fsamp, thisRfs, 'nearest')
                else:
                    # H-pol propagation factor composed of:
                    #    1. Interference region
                    F = interp1(rngSamp, FinterfSamp, thisRfs, 'nearest')
            case 'V':
                # V-pol propagation factor composed of:
                #    1. Interference region
    
                # Calculate direct ray elevation angles
                thetad, _, _, minHt, flipGeom = directRayElAng(thisRfs, anht, tgtht, Rhoriz, Re)
    
                # Shift antenna pattern
                anpat, patel, tiltOut = antennaPatternShift(anpat, patel, tilt, flipGeom)
    
                # Interference propagation factor 
                F, _, _ = propfactorinterf(freq, minHt, thetad, epsc, pol, anpat, patel, tiltOut, htsd, beta0, vegTypeP, Re)
            case _:
                raise ValueError("Polarization should be either 'H' or 'V'")
    
        # Set NaN values to a small number (NaNs are generally at ranges > Rhoriz)
        FdB[idx] = mag2db(F)
        FdB[isnan(FdB)] = -inf
    
    return FdB

def effearthradius(
    refgrad = -39e-9,
    R = None,
    ha = None,
    ht = None,
    ns = 313,
    altbp = None,
    npb = None,
):
    '''
    R: line-of-sight range to target
    ha: radar altitude above mean sea level
    ht: target altitude above mean sea level
    ns: "SurfaceRefractivity"
    altbp: "BreakPointAltitude"
    npb: "BreakPointRefractivity"
    
    %effearthradius    Effective earth radius
    %   Re = effearthradius returns the effective radius of a spherical earth
    %   (in meters). The effective radius is calculated using a refractivity
    %   gradient of -39e-9, which results in approximately 4/3 of the real
    %   earth radius.
    %
    %   Re = effearthradius(RGradient) specifies the refractivity gradient
    %   RGradient (in units of 1/meter).
    %
    %   Re = effearthradius(R,HA,HT) returns the effective Earth radius, Re, in
    %   meters using the average radius of curvature method. The input R is
    %   either a scalar or the M-length line-of-sight range to the target in
    %   meters, HA is either a scalar or the M-length radar mean sea level
    %   (MSL) altitude in meters, HT is either a scalar or the M-length target
    %   MSL altitude in meters.
    %
    %   Re = effearthradius(R,HA,HT,'SurfaceRefractivity',NS) calculates the
    %   effective Earth radius using the provided scalar surface refractivity
    %   NS in N-units. Defaults to 313.
    %
    %   Re = effearthradius(R,HA,HT,...,'BreakPointAltitude',ALTBP) defines the
    %   altitude of the convergence point as a scalar in meters. Defaults to
    %   12192 meters (40 kft) when any of the input altitudes are greater than
    %   9144 meters (30 kft). Otherwise defaults to 9144 meters (30 kft). The
    %   convergence point improves errors at higher altitudes by sacrificing
    %   fidelity at lower altitudes.
    %
    %   Re = effearthradius(R,HA,HT,...,'BreakPointRefractivity',NBP) defines
    %   the refractivity of the convergence point as a scalar in N-units.
    %   Defaults to 66.65 N-units when any of the input altitudes are greater
    %   than 9144 meters (30 kft). Otherwise defaults to 102.9 N-units. The
    %   convergence point improves errors at higher altitudes by sacrificing
    %   fidelity at lower altitudes.
    %
    %   [Re,k] = effearthradius(...) also outputs the effective radius factor k
    %   with the default approximately equal to 4/3.
    %
    %   The effective Earth radius factor k is an approximation for refraction
    %   effects in the troposphere. It ignores other types of propagation
    %   phenomena such as ducting.
    %
    %   Often, refraction is accounted for by assuming an effective radius
    %   factor equal to 4/3. However, at long ranges and with shallow angles,
    %   the k factor can deviate greatly from the 4/3 approximation. No
    %   atmospheric refraction results in a value of k = 1. A value of k =
    %   infinity is a flat Earth.
    %
    %   % Example 1:
    %   %   Calculate the effective earth radius.
    %   re = effearthradius
    %
    %   % Example 2:
    %   %   Calculate the effective Earth radius Re for a radar at 0 meters and
    %   %   a target at 7.62e3 meters (25 kft) at a range of 250 km.
    %   R = 250e3 % Line-of-sight range (m) 
    %   ha = 0 % Radar altitude (m) 
    %   ht = 7.62e3 % Target altitude (m)
    %   re = effearthradius(R,ha,ht)
    % 
    %   % Example 3:
    %   %   Compare the average curvature method to the 4/3 earth model. 
    %   R = 100e3 % Line-of-sight range (m) 
    %   ha = (0.25:0.01:15).*1e3 % Radar altitude (m) 
    %   ht = 0 % Target altitude (m)
    %   [~,kAvgCurv] = effearthradius(R,ha,ht)
    % 
    %   % Plot results
    %   semilogy(kAvgCurv,ha*3.2808*1e-3)
    %   hold on
    %   plot(4/3*ones(1,2),[10^-3 10^3],'--k')
    %   grid on
    %   axis([1.05 1.4 10^0 10^2])
    %   legend('Average Curvature','4/3 Earth')
    %   xlabel('Effective Earth Radius Factor k')
    %   ylabel('Altitude (kft)')
    %
    %   See also horizonrange, depressionang, grazingang.
    
    %   Copyright 2010-2020 The MathWorks, Inc.
    
    %   Reference
    %   [1] Doerry, A. W. "Earth Curvature and Atmospheric Refraction Effects
    %       on Radar Signal Propagation." Sandia National Laboratories,
    %       SAND2012-10690, Jan. 2013.
    %   [2] Long, Radar Reflectivity of Land and Sea, 3rd Ed. Artech House, 
    %       2001.
    %   [3] Mahafza, Radar Signal Analysis and Processing Using MATLAB, CRC 
    %       Press, 2009.
    %   [4] Skolnik, Introduction to Radar Systems, 3rd Ed. McGraw Hill, 2002.
    %   [5] Ward, Space-Time Adaptive Processing for Airborne Radar, 
    %       Lincoln Lab Technical Report, 1994.
    
    %   In the original syntax, with no refraction, the curvature is 1/R, which
    %   means that the horizontal beam is always parallel to the earth.
    %   However, with the change of refractivity for each layer, it's no longer
    %   true. The effective earth radius is the radius that the modified
    %   curvature is 1/Re.
    '''
    if (ha is not None) and (ht is not None):
        if altbp is None:
            if any([all(x > 9144 for x in ht), all(x > 9144 for x in ha)]): # Interested in range of altitudes between 0 and 50 kft (15240 m)
                altbp = 12192 # hb, m
            else: # Interested in range of altitudes between 0 and 30 kft (9144 m)
                altbp = 9144 # hb, m
        if npb is None:
            if any([all(x > 9144 for x in ht), all(x > 9144 for x in ha)]): # Interested in range of altitudes between 0 and 50 kft (15240 m)
                npb = 66.65 # Nb, N-units
            else: # Interested in range of altitudes between 0 and 30 kft (9144 m)
                npb = 102.9 # Nb, N-units

    # Calculate Earth radius and k factor
    if any(x is None for x in [R, ha, ht, altbp, npb]):
        k = 1 / (1 + Rearth * refgrad)
        Re = k * Rearth
    else:
        Re, k = effearthmethods(R, ha, ht, ns, altbp, npb)
    
    return Re, k

def effearthmethods(R, ha, ht, Ns, hb, Nb):
    '''
    %   This class is for internal use only. It may be removed in the future.
    
    % effearthmethods    Calculate effective earth radius
    %   [Re,k] = effearthmethods(R,HA,HT,METHOD,Ns,hb,Nb) returns the effective
    %   Earth radius in m using the Average Curvature method. 
    %
    %   The inputs are as follows:
    %      - R is either a scalar or the M-length line-of-sight range to the 
    %        target in m
    %      - HA is either a scalar or the M-length radar mean sea level
    %        (MSL) altitude in m
    %      - HT is either a scalar or the M-length target MSL altitude in m
    %      - Ns defines scalar surface refractivity in N-units
    %      - hb defines the altitude of the convergence point as a scalar in m
    %      - Nb defines the refractivity of the convergence point as a scalar 
    %        in N-units
    %
    %   This function returns Re, the effective radius of a spherical earth in
    %   meters and k, the effective radius factor. 
    %
    %   Copyright 2019-2020 The MathWorks, Inc.
    
    %   Reference 
    %   [1] Doerry, A. W. "Earth Curvature and Atmospheric Refraction Effects
    %       on Radar Signal Propagation." Sandia National Laboratories,
    %       SAND2012-10690, Jan. 2013.
    '''    
    # Define iterative calculation stop criteria
    tol = 1e-6
    maxIterations = 10
    
    # Calculate break point
    numHa = numel(ha)
    numHt = numel(ht)
    if numHa == numHt:
        heights = np.array([[ha], [ht]])
    else:
        if numHt == 1:
            heights = np.vstack((ha, ht * ones((1, numHa))))
        else:
            heights = np.vstack((ha * ones((1, numHt)), ht))
            
    hmin, _ = matlabmin(heights, [], 1)
    hmax, _ = matlabmax(heights, [], 1)
    dh = hmax - hmin
    Hb = (hb - hmin) / log(Ns / Nb)
    
    # Estimate initial grazing angle
    k0 = 4 / 3
    gAng0 = asind((dh) / R * (1 + dh / (2 * (k0 * Rearth + hmin))) - R / (2 * (k0 * Rearth + hmin)))
    N = size(gAng0, 2)
    
    # Initialize iterative search
    nIterations = 1
    gAngIt = zeros((2, N))
    gAngIt[1, :] = gAng0
    kIt = zeros((2, N))
    kIt[1, :] = k0
    ctol = Inf(1, N)
    
    # Iterate to determine k and grazing angle
    while np.any(x > tol for x in ctol) and nIterations <= maxIterations:
        # Update k
        b = exp(dh / Hb)
        c = 10 ** -6 * Ns * cosd(gAngIt[1, :]) * Rearth / Hb
        
        # Use "Average Radius of Curvature" approximation method
        d = dh / Hb / (b - 1)
        kIt[0, :] = 1 / (1 - c * d) # Eqn 37
        
        # Update grazing angle
        gAngIt[0, :] = asind((dh) / R * (1 + dh / (2 * (kIt[0, :] * Rearth + hmin))) - R / (2 * (kIt[0, :] * Rearth + hmin)))
        
        # Calculate tolerance
        ctol = abs(kIt[0, :] - kIt[1, :])
        
        # Prepare for next iteration
        gAngIt[1, :] = gAngIt[0, :]
        kIt[1, :] = kIt[0, :]
        nIterations = nIterations + 1
    
    # Define outputs
    k = kIt[1, :]
    Re = k * Rearth

    return Re, k

def radarpropfactor_getRange(rfs, it, numTgt):
    if numTgt == 1:
        idx = list(range(0, np.size(rfs)))
    else:
        idx = it

    thisRfs = rfs[idx]

    return thisRfs, idx

def radarpropfactor_getRangeSamples(rfs, rmax, rmin):
    if (np.size(rmax) == 1) and isinstance(rmax, np.ndarray):
        rmax = rmax[0]
    if (np.size(rmin) == 1) and isinstance(rmin, np.ndarray):
        rmin = rmin[0]
    
    # Calculate range samples to standardize results    
    exp0 = matlabmax(floor(log10(rmin)), 0)
    expe = ceil(log10(rmax)).astype(int)
    
    # Returns range sample points in range [0 rmax]
    numPointsPerDecade = 1e4
    numDecades = ceil(np.absolute(expe - exp0) / 10)
    
    rngSamp = logspace(exp0, expe, np.array(numPointsPerDecade * numDecades).astype(int))
    
    # Calculation ranges = sorted unique values including samples defined above
    # and user requested values
    rngSampOut, _ = unique([*list(rngSamp[rngSamp > rmin]), *list(rfs)])
    
    return rngSampOut

def directRayElAng(R, anht, tgtht, Rhoriz, Re):
    # Calculate direct ray elevation angles

    if isscalar(R):
        R = [R]

    # Limit calculations for interference region to horizon
    # idxRi = [i for i, x in enumerate(R) if x < Rhoriz]
    idxRi = np.where(R <= Rhoriz)[0] # Do not perform interference calculations beyond horizon
    Ri = np.full(np.shape(R), nan)
    Ri[idxRi] = R[idxRi]
    
    # Do not permit calculations where range is less than minimum range
    maxHt = max(tgtht, anht)
    minHt = min(tgtht, anht)
    val = (2 * Re * (maxHt - minHt) + maxHt ** 2 - minHt ** 2 - Ri ** 2) / (2 * (Re + minHt) * Ri)
    idxGt = abs(val) > 1
    val[idxGt] = nan # Limit the arcsine value to 1
    
    # The original Blake paper was written for positive elevation angles.
    # radarpropfactor handles negative elevation angles by reversing the
    # geometry. Change the sign of the tilt angle in the case of a negative
    # elevation geometry (i.e., the radar is higher than the target). 
    flipGeom = False
    if anht > tgtht:
        # Geometry is flipped
        flipGeom = True
    
    # Calculate angle from sensor to target
    thetad = arcsin(val) # Ref 2, Eqn 8.47 (rad)
    return thetad, Ri, maxHt, minHt, flipGeom

def antennaPatternShift(anpatIn, patelIn, tiltIn, flipGeom):
    if flipGeom:
        # Geometry is flipped. Shift antenna pattern for negative elevation
        # geometries.
        patel = -(patelIn + tiltIn)
        patel = mod(patel + 90, 180) - 90 # Wrap to +- 90
        patel, idxUnique = unique(np.ravel(patel))
        patel = np.sort(patel)
        anpat = anpatIn[idxUnique]
        tiltOut = 0
    else:
        anpat = anpatIn
        patel = patelIn
        tiltOut = tiltIn 

    return anpat, patel, tiltOut

def propfactorinterf(freq, minHt, thetad, epsc, pol, anpatIn, patelIn, tilt, htsd, beta0, vegType, Re):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %propfactorinterf Radar propagation factor for the interference region
    %   [FI,RI,PD] = propfactorinterf(FREQ, MINHT, THETAD, EPSC, POL, ANPAT,
    %   PATEL, TILT, HTSD, BETA0, VEGTYPE, RE) calculates the radar propagation
    %   factor in the interference region in linear units.
    %
    %   The inputs are as follows:
    %      FREQ    = Scalar input frequency (Hz)
    %      MINHT   = Minimum height of either radar or target (m)
    %      THETAD  = Direct ray elevation angles (rad)
    %      EPSC    = Scalar complex permittivity
    %      POL     = Polarization 'H'|'V'
    %      ANPAT   = Antenna pattern (linear units)
    %      PATEL   = Antenna pattern elevation angles (deg) 
    %      TILT    = Tilt angle (deg) 
    %      HTSD    = Surface height standard deviation (m) 
    %      BETA0   = Surface slope (deg) 
    %      VEGTYPE = Vegetation type 'Trees'|'Weeds'|'Brush'|'Grass'|'None'
    %      RE      = Effective Earth radius (m) 
    %   See the calling function for additional information. 
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also radarpropfactor.
    
    %   Copyright 2020-2022 The MathWorks, Inc.  
    
    %   References
    %   [1] Blake, L.V. "Machine Plotting of Radar Vertical-Plane Coverage
    %       Diagrams." Naval Research Laboratory, 1970 (NRL Report 7098).
    %   [2] Barton, David K. Radar Equations for Modern Radar. Norwood, MA:
    %       Artech House, 2013.
    '''
    
    # Calculate wavelength
    wv, _ = freq2wavelen(freq) # Wavelength (m)
    
    # Set angles 
    anpat = singrowvec(anpatIn)
    patel = deg2rad(singrowvec(patelIn)) # Antenna pattern angles (rad)
    tilt = deg2rad(tilt)          # Tilt of antenna (rad)
    beta0 = deg2rad(beta0)        # Surface slope (rad)
    phi = sqrt((tan(thetad) / 3) ** 2 + 2 * minHt / (3 * Re)) - tan(thetad) / 3 # Phase angle of complex reflection coefficient (rad) Eq 25, ref 1
    
    # Compute path difference
    pdCE = pathdiffCE(thetad, minHt, phi, Re) # Path difference for curved earth, Eqn 24, ref 1
    pdFE = pathdiffFE(thetad, minHt)        # Path difference for flat earth, Eqn 23, ref 1
    idxH = np.absolute(pdCE - pdFE) < 0.01 * wv  # Higher elevation, see ref 1, subroutine LOBES
    idxL = ~idxH                           # Lower elevation, subroutine LOBES
    pd = zeros(np.size(pdCE))
    pd[idxL] = pdCE[idxL]
    pd[idxH] = pdFE[idxH]
    graz = thetad + phi                    # Grazing angle (rad) 
    thetar = -graz                         # Angle of reflected wave, Fig 1 (rad) 
    thetar[idxH] = -thetad[idxH]
    
    # Calculate reflection coefficient
    rhoref, phiref = reflectioncoeff(epsc, graz, pol) # Eqn 29, 30, ref 1
    
    # Determine vegetation factor
    rhov = vegetationfactor(freq, graz, vegType) # Vegetation factor, ref 2
    
    # Determine effective height
    # Equation 2 in command-line help: Eqn 20, ref 1
    # Assumes
    #    htsd = Standard deviation of height = RMS
    #         = (crest-to-trough height)/(2*sqrt(2))
    #         = (amplitude)/(sqrt(2))
    effHtScale = np.absolute(2 * graz / beta0)
    if any(effHtScale < 1): # Condition for shadowing
        effectiveHt =  htsd * (effHtScale) ** (0.2) # Eqn 8.59, ref 2
    else:
        effectiveHt = htsd * ones(np.size(graz))
    
    # Calculate specular scattering (i.e., roughness)
    r = roughness(wv, effectiveHt, graz) # Eqn 20, ref 1
    
    # Calculate divergence factor
    D = divfactor(minHt, thetad, Re) # Eqn 27, ref 1
    D[idxH] = 1
    
    # Calculate phase difference, Eqn 9, ref 1
    ratio = pd / wv + phiref / (2 * pi) # ref 1, frplot
    ratio1 = rem(ratio, 1)
    alpha = 2 * pi * ratio1
    
    # Propagation factor Eqn 8, ref 1
    theta1 = mod(((thetad - tilt) + pi / 2), pi) - pi / 2 # Wrap to +- 90
    theta2 = mod(((thetar - tilt) + pi / 2), pi) - pi / 2 # Wrap to +- 90
    apatd = antfactor(anpat, patel, theta1)
    apatr = antfactor(anpat, patel, theta2)
    grc = (r * D * rhoref * rhov * apatr / apatd) # Eqn 7, ref 1
    Fi = np.absolute(apatd) * (sqrt(1 + grc ** 2 + 2 * grc * cos(alpha))) # Eqn 8, ref 1
    
    return Fi, pd, wv
    
def antfactor(anpat, patel, theta):
    # antresp Returns interpolated antenna response at theta based on
    # provided response samples anpat and patangles
    antresp = interp1(singrowvec(patel), singrowvec(anpat), singrowvec(theta))
    
    return antresp
    
def pathdiffFE(thetad, ha):
    # Path difference for flat Earth assumption
    # Ref: NRL Report 7098
    delFE = 2 * ha * sin(thetad)
    
    return delFE
    
def pathdiffCE(thetad, ha, phi, Re):
    # Path difference for curved Earth assumption
    # Ref: NRL Report 7098
    delCE = sqrt(ha ** 2 + Re * (Re + ha) * phi ** 2) * sin(thetad + phi) ** 2 * 2
    
    return delCE
    
def roughness(wv, htsd, graz):
    # roughness Returns roughness factor of surface assuming Gaussian height
    # distribution. 
    # wv Wavelength (m)
    # HTSD Standard deviation of height (m)
    # GRAZ Grazing angle (rad) 
    # Ref NRL Report 7098 eq 20
    r = exp(-2 * ((2 * pi * htsd * sin(graz)) / wv) ** 2)
    return r

def reflectioncoeff(epsc, graz, pol):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %reflectioncoeff Calculate reflection coefficient
    %   [RHOREF, PHIREF] = reflectioncoeff(EPSC,GRAZ,POL) calculates the
    %   reflection coefficient. RHOREF is the magnitude of the complex
    %   reflection coefficient and PHIREF is the phase of the complex
    %   reflection coefficient. 
    %
    %   EPSC is expected to be the scalar complex relative permittivity, GRAZ
    %   is a vector of grazing angles in radians, and POL is either 'H'|'V' for
    %   horizontal and vertical polarization, respectively. 
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also radarpropfactor.
    
    %   Copyright 2020 The MathWorks, Inc.
    
    %   References
    %   [1] Blake, L.V. "Machine Plotting of Radar Vertical-Plane Coverage
    %       Diagrams." Naval Research Laboratory, 1970 (NRL Report 7098).
    %   [2] Barton, David K. Radar Equations for Modern Radar. 1st edition.
    %       Norwood, MA: Artech House, 2013.
    '''    
    match pol:
        case 'H':
            # Horizontal polarization
            refc = (sin(graz) - sqrt(epsc - cos(graz) ** 2)) / (sin(graz) + sqrt(epsc - cos(graz) ** 2))

        case 'V':
            # Vertical polarization
            refc = (epsc * sin(graz) - sqrt(epsc - cos(graz) ** 2)) / (epsc * sin(graz) + sqrt(epsc - cos(graz) ** 2))

        case _:
            raise ValueError("Polarization needs to be either 'H' or 'V'")

    rhoref = abs(refc)
    phiref = -angle(refc)
    
    return rhoref, phiref

def seaComplexPermittivity(freq):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %seaComplexPermittivity Blake's complex permittivity model of sea water
    %   EPSC = seaComplexPermittivity(FREQ) calculates the scalar complex
    %   permittivity EPSC based on the scalar input frequency FREQ in Hz. 
    %
    %   Model has a valid range of 100 MHz to 10 GHz.
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also radarpropfactor.
    
    %   Copyright 2020 The MathWorks, Inc.
    
    %   References
    %   [1] Blake, L.V. "Machine Plotting of Radar Vertical-Plane Coverage
    %       Diagrams." Naval Research Laboratory, 1970 (NRL Report 7098).
    '''    
    # NRL Report 7098 Eqn 32
    freqMhz = np.array(freq / 1e6) # Convert to MHZ, used formula is for MHz
    
    # Initialize
    epsr = zeros(np.shape(freq))
    sigma = zeros(np.shape(freq))
    
    # Range 1: Frequency <= 1.5 GHz
    idx = freqMhz <= 1500
    if idx:
        if np.size(idx) > 1:
            epsr[idx] = 80
            sigma[idx] = 4.3
        else:
            epsr = 80
            sigma = 4.3
    
    # Range 2: Frequency > 1.5 GHz and <= 3 GHz
    idx = ((freqMhz > 1500) & (freqMhz <= 3000))
    if idx:
        if np.size(idx) > 1:
            epsr[idx] = 80 - 0.00733 * (freqMhz[idx] - 1500)
            sigma[idx] = 4.30 + 0.00148 * (freqMhz[idx] - 1500)
        else:
            epsr = 80 - 0.00733 * (freqMhz - 1500)
            sigma = 4.30 + 0.00148 * (freqMhz - 1500)
    
    # Range 3: Frequency > 3 GHz (valid to 10 GHz) 
    idx = ((freqMhz > 3000) & (freqMhz <= 10000))
    if idx:
        if np.size(idx) > 1:
            epsr[idx] = 69 - 0.00243 * (freqMhz[idx] - 3000)
            sigma[idx] = 6.52 + 0.001314 * (freqMhz[idx] - 3000)
        else:
            epsr = 69 - 0.00243 * (freqMhz - 3000)
            sigma = 6.52 + 0.001314 * (freqMhz - 3000)
    
    # Calculate wavelength
    wv, _ = freq2wavelen(freq)
    
    # Calculate permittivity
    epsc = epsr - 1j * 60 * np.array(wv) * np.array(sigma)
    
    return epsc
    
def searoughness(n, scaletype = "SeaState"):
    '''
    %searoughness Surface height standard deviation for sea
    %   HGTSD = searoughness(N) outputs the standard deviation of surface
    %   height for the specified sea state number as a scalar in m (i.e.,
    %   surface roughness).
    %
    %   The sea state N is a nonnegative scalar integer. Acceptable sea state
    %   values are in the range from [0 - 8], inclusive.
    %
    %   HGTSD = searoughness(...,'ScaleType',SCALETYPE) specifies the
    %   scale type as either 'SeaState' | 'WindScale'. The sea state model
    %   accepts an input N that is a nonnegative scalar integer in the range [0
    %   - 8]. The wind scale is the Beaufort Wind Scale model whose inputs are
    %   positive scalar integers in the range of [1 - 9]. Defaults to
    %   'SeaState'.
    %
    %   [HGTSD,BETA0,VW] = searoughness(...) additionally outputs the
    %   scalar slope of the sea type BETA0 in degrees and the scalar wind
    %   velocity VW in m/s. Note that BETA0 is 1.4 times the RMS surface slope.
    %
    %   % Example:
    %   %   Obtain the surface height standard deviation assuming a sea state
    %   %   of 2.
    %   hgtsd = searoughness(2)
    %
    %   See also landroughness, seareflectivity, landreflectivity,
    %   clutterSurfaceRCS, radarpropfactor, radarvcd, blakechart.
    
    %   Copyright 2020 The MathWorks, Inc.
    
    %   Reference
    %   [1] Barton, David K. Radar Equations for Modern Radar. 1st edition.
    %       Norwood, MA: Artech House, 2013.
    '''    
    # Validate sea type
    if scaletype == "SeaState": # Sea State must be 0-8
        ss = n   # Sea state
    elif scaletype == "WindScale": # Beaufort Wind Scale must be 1-9
        ss = n - 1 # Equivalent sea state
    else:
        raise ValueError("scaletype must be either 'SeaState' or 'WindScale'")
    
    # Setup table (Ref 1, Table 9.1)
    windVel = [1.5, 2.6, 4.6, 6.7, 8.2, 10.8, 13.9, 19.0, 28.8] # m/s
    rmsHt = [0.01, 0.03, 0.10, 0.24, 0.38, 0.57, 0.91, 1.65, 2.50] # m
    slope = [0.055, 0.063, 0.073, 0.080, 0.085, 0.091, 0.097, 0.104, 0.116] # rad
    
    # Set Gaussian model outputs
    idx   = ss              # Sea state starts at 0
    vw    = windVel[idx]        # m/s
    hgtsd  = rmsHt[idx]         # m
    beta0 = np.degrees(slope[idx]) # deg

    return hgtsd, beta0, vw

def sincpattern(hpbw, ang = np.linspace(-90, 90, 361), K = 1.39157):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %sincpattern Theoretical sinc antenna pattern
    %   APAT = sincpattern(HPBW,ANG,K) calculates the theoretical sinc antenna
    %   pattern, which is of the form sin(x)/x. 
    % 
    %   HPBW = Half-power beamwidth in azimuth or elevation (degrees)
    %   ANG  = Azimuth or elevation angles (degrees)
    %   K    = Beam factor (typically set to 1.39157 so that APAT is = 0.7071 =
    %          1/sqrt(2) when ANG = HPBW/2 (unitless))
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also radarpropfactor, radarvcd, phased.SincAntennaElement.
    
    %   Copyright 2021 The MathWorks, Inc.
    
    %   References
    %   [1] Blake, L.V. "Machine Plotting of Radar Vertical-Plane Coverage
    %       Diagrams." Naval Research Laboratory, 1970 (NRL Report 7098).
    '''
    
    # Sinc antenna pattern
    apat = sinc(K / np.degrees(sin(np.radians(np.array(hpbw) / 2))) * np.degrees(sin(np.radians(ang))) / pi) # Eqns 38 - 40 
    
    return apat, ang

def vegetationfactor(freq, graz, vegType):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %vegetationFactor    Calculate the absorption factor for vegetation 
    %   RHOV = vegetationFactor(FREQ,GRAZ,VEGTYPE) calculates the
    %   absorption factor for vegetation. FREQ is a J-length vector of
    %   frequencies in Hz, GRAZ is a K-length vector of grazing angles in
    %   radians, and VEGTYPE is one of either 'Trees' | 'Weeds' | 'Brush' |
    %   'Grass' | 'None'.
    %
    %   The output RHOV is a JxK matrix of vegetation factors where J
    %   corresponds to the number of frequencies and K corresponds to the
    %   number of grazing angles reported in linear units. 
    %
    %   In the case of VEGTYPE set to 'Trees' | 'Weeds' | 'Brush', there is an
    %   assumption of dense vegetation. The case of VEGTYPE set to 'Grass'
    %   assumes thin grass. 
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also radarpropfactor. 
    
    %   Copyright 2020 The MathWorks, Inc.
    
    %   Reference
    %   [1] Barton, David K. Radar Equations for Modern Radar. 1st edition.
    %       Norwood, MA: Artech House, 2013.
    '''
    
    # Calculate wavelength (m)
    wv, _ = freq2wavelen(freq)
    
    # Coefficients a and b (See section 8.3.3, pg 295 and Supplementary
    # Material sheet 9-1)
    match vegType:
        case 'Trees':
            a = 0.0316
            b = 5
        case 'Brush' | 'Weeds':
            a = 0.316
            b = 3
        case 'Grass':
            a = 3.16
            b = 1
        case _: # 'None'
            a = 0 
            b = 0 
    
    # Calculate vegetation factor 
    expVal = (-b * sin(graz)) / wv
    rhov = ((1 - sqrt(a * wv)) * exp(expVal)) + sqrt(a * wv) # Eqn 8.61 
    rhov[rhov < 0] = 0
    rhov[rhov > 1] = 1 

    return rhov

def divfactor(ha, thetad, Re):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %divfactor Returns divergence factor 
    %   D = divfactor(HA,THETAD,RE) returns the divergence factor in linear
    %   units. HA is the scalar antenna height in m, THETAD is the elevation
    %   angle of the direct ray path in radians, and RE is the scalar effective
    %   Earth radius in m.
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also radarpropfactor.
    
    %   Copyright 2020-2021 The MathWorks, Inc.
    
    %   References
    %   [1] Blake, L.V. "Machine Plotting of Radar Vertical-Plane Coverage
    %       Diagrams." Naval Research Laboratory, 1970 (NRL Report 7098).
    
    '''
    
    # Eqn 27, 28
    if ha < sqrt(eps()):
        D = ones(np.size(thetad))
    else:
        zeta = sqrt(Re / (2 * ha)) * tan(thetad)
        idxNaN = isinf(zeta)
        val = ones(np.size(thetad))
        val[~idxNaN] = (1 + (2 * zeta[~idxNaN] / sqrt(zeta[~idxNaN] ** 2 + 3))) / 3
        D = sqrt(val)

    return D

def propfactordiff(Rin, freq, minHt, maxHt, Re, n0):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %propfactordiff Radar propagation factor for the diffraction region
    %   FD0 = propfactordiff(R,FREQ,MINHT,MAXHT,RE,NO) calculates the radar
    %   propagation factor in the diffraction region in linear units.
    %
    %   R is a row vector of ranges in m, FREQ is the scalar input frequency in
    %   Hz, MINHT is the minimum altitude in the geometry in m, MAXHT is the
    %   maximum altitude in the geometry in m, RE is the scalar effective Earth
    %   radius in m, and N0 is the scalar surface refractive index.
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also radarpropfactor.
    
    %   Copyright 2020 The MathWorks, Inc.
    
    %   References
    %   [1] Barton, David K. Radar Equations for Modern Radar. 1st edition.
    %       Norwood, MA: Artech House, 2013.
    '''
    
    # Convert to row vector
    R = singrowvec(Rin) 
    
    # Wavelength
    wv, _ = freq2wavelen(freq)  
    
    # Natural units are given as a function of wavelength and the effective
    # Earth radius
    L = ((Re ** 2) * wv / (pi * n0)) ** (0.333)   # m
    H = (Re * wv ** 2 / (8 * pi ** 2 * n0)) ** (0.333) # m
    
    # Supporting variables for Eqn 8.63
    X = R / L      # Range in natural units
    Zr = minHt / H # Terminal height sensor in natural units
    Zt = maxHt / H # Terminal height target in natural units
    
    # Attenuation factor
    V = rngfactor(X)
    
    # Calculate height-gain factors 
    Ur = htgainfactor(Zr) # Height-gain factor for sensor
    Ut = htgainfactor(Zt) # Height-gain factor for target 
    
    # Eqn 8.63
    #    Fd0 = 2*sqrt(pi*X)*exp(-2.02*X)*|U(Zr)*U(Zt)|
    #        = V(X)*|U(Zr)*U(Zt)|
    # To convert back to linear units use 10^(Fd0/20).
    Fd0 = V + Ur + Ut # dB
    Fd0 = db2mag(Fd0)  # Linear units
    
    # Eliminate invalid high values that result from use of a single mode 
    Fd0 = sqrt(Fd0 ** 2 / (Fd0 ** 2 + 1)) # Eqn 8.69 (linear units)
    
    return Fd0
    
def htgainfactor(Z):
    # Set calculation types
    idxLt = np.where(Z <= 0.6)
    idxBtwn = np.where((Z < 1) & (Z > 0.6))
    idxGt = np.where(Z >= 1)
    
    # Define U (Eqn 8.66)
    # Produces U in dB. To convert back to linear units use 10^(U/20).
    U = zeros(np.size(Z))
    
    if np.size(idxLt) > 0:
        if np.size(idxLt) > 1:
            U[idxLt] = 20 * log10(Z[idxLt])
        else:
            U = 20 * log10(Z)
    if np.size(idxBtwn) > 0:
        if np.size(idxBtwn) > 1:
            U[idxBtwn] = -4.3 + 51.04 * (log10(Z[idxBtwn] / 0.6)) ** (1.4)
        else:
            U = -4.3 + 51.04 * (log10(Z / 0.6)) ** (1.4)
    if np.size(idxGt) > 0:
        if np.size(idxGt) > 1:
            U[idxGt] = 19.84728 * (Z[idxGt] ** 0.47 - 0.9)
        else:
            U = 19.84728 * (Z ** 0.47 - 0.9)

    return U
    
def rngfactor(X):
    # Define V (Eqn 8.67)
    # Produces V in dB. To convert back to linear units use 10^(V/20).
    V = 10.99209864 + 10 * log10(X) - 17.545497 * X 

    return V

def propfactorall(Rin, Finterfin, Fdiffin, Rdelta, Rhoriz):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %propfactorall Composite radar propagation factor for all regions
    %   F = propfactorall(R,FINTERF,FDIFF,RDELTA,RHORIZ) calculates
    %   the radar propagation factor in the intermediate region, an
    %   interpolation between the interference and diffraction regions. Outputs
    %   the composite propagation factor composed of the interference,
    %   intermediate, and diffraction regions in linear units.
    %
    %   The inputs are as follows: 
    %      R      = Row vector of ranges (m)
    %      Fi     = Interference region propagation factor (linear units)
    %      Fd     = Diffraction region propagation factor (linear units)
    %      Rdelta = Delta range, beginning of intermediate region (m)
    %      Rhoriz = Horizon range, beginning of diffraction region (m) 
    %
    %   This function assumes the propagation factor involves the following 3
    %   regions: 
    %      1. Interference region: <= Rdelta
    %      2. Intermediate region: > Rdelta and < Rhoriz
    %      3. Diffraction region: >= Rhoriz
    % 
    %   See the calling function for additional information. 
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also radarpropfactor.
    
    %   Copyright 2020 The MathWorks, Inc.
    
    %   References
    %   [1] Blake, L.V. "Machine Plotting of Radar Vertical-Plane Coverage
    %       Diagrams." Naval Research Laboratory, 1970 (NRL Report 7098).
    %   [2] Barton, David K. Radar Equations for Modern Radar. 1st edition.
    %       Norwood, MA: Artech House, 2013.
    '''
        
    # Convert to row vectors
    R = singrowvec(Rin)
    Finterf = singrowvec(Finterfin)
    Fdiff = singrowvec(Fdiffin)
    
    # Convert propagation factors to dB
    FddB = mag2db(Fdiff)
    FidB = mag2db(Finterf)
    
    # Initialize output
    FdB = Nan(np.size(R))
    
    # Determine indices for interpolation
    _, idxHoriz = matlabmin(np.absolute(R - Rhoriz), [], 2)
    _, idxDelta = matlabmin(np.absolute(R - Rdelta), [], 2)
    idxII = list(range(idxDelta, (idxDelta + 2)))
    idxID = list(range((idxHoriz - 1), idxHoriz + 1))
    
    # Check if less than 4 unique points exist
    errorCond = numel(unique([idxII, idxID])[0]) != 4 # Need 4 unique points 
    if errorCond:
        # Do not interpolate between regions. Exit function. 
        F = Finterf
        return
    
    # Interpolate intermediate region, reference 2, section 8.6.1
    rngInterp = np.hstack((R[idxII], R[idxID]))
    Finterp = np.hstack((FidB[idxII], FddB[idxID]))
    idxInt = list(range((idxDelta + 1), idxHoriz))
    FdB[0:idxDelta] = FidB[0:idxDelta] # Interference region, Eqn 8.78, ref 2
    FdB[idxInt] = interp1(rngInterp, Finterp, R[idxInt], 'pchip') # 4 point, cubic interpolation, ref 1, frplot
    FdB[idxHoriz:] = FddB[idxHoriz:] # Diffraction region, Eqn 8.78, ref 2
    
    # Convert to linear units
    F = db2mag(FdB)
    
    return F

def effbeamwidth(thetaTx, thetaRx):
    '''
    %effbeamwidth Effective beamwidth
    %   THETAEFF = effbeamwidth(THETATX,THETARX) calculates the positive
    %   effective beamwidth THETAEFF in degrees.
    %
    %   The inputs THETATX and THETARX are either scalar or length-J vector
    %   positive azimuth or elevation beamwidths, where Tx represents the
    %   transmitter beamwidth (degrees) and Rx represents the
    %   receiver beamwidth (degrees). The effective beamwidth can be
    %   used for the calculation of the clutter illumination area or the
    %   effective angular resolution in a simulated radar. Note that the units
    %   of THETATX and THETARX must be the same.
    %
    %   % Example 1:
    %   %   Calculate the effective azimuth beamwidth in degrees for a radar 
    %   %   system that has a transmitter azimuthal beamwidth of 60 degrees 
    %   %   and a receiver azimuthal beamwidth of 5 degrees. 
    %   thetaAzTx = 60 % Transmitter azimuthal beamwidth (deg) 
    %   thetaAzRx = 5  % Receiver azimuthal beamwidth (deg)
    %   thetaAzEff = effbeamwidth(thetaAzTx,thetaAzRx)
    %
    %   % Example 2:
    %   %   Calculate the effective elevation beamwidths in degrees for radar 
    %   %   systems that have transmitter elevation beamwidths of 10 and 15  
    %   %   degrees, respectively, and receiver elevation beamwidths of 2 
    %   %   and 5 degrees, respectively. 
    %   thetaElTx = [10 15] % Transmitter elevation beamwidths (deg) 
    %   thetaElRx = [2 5]  % Receiver elevation beamwidths (deg)
    %   thetaElEff = effbeamwidth(thetaElTx,thetaElRx)
    %
    %   See also ap2beamwidth, beamwidth2ap, beamwidth2gain. 
    
    %   Copyright 2020-2024 The MathWorks, Inc.
    
    %   Reference
    %   [1] Barton, David K. Radar Equations for Modern Radar. 1st edition.
    %       Norwood, MA: Artech House, 2013.
    '''
    
    numThetaTx = numel(thetaTx)
    numThetaRx = numel(thetaRx)

    if not isinstance(thetaTx, np.ndarray):
        if (numThetaTx > 1):
            thetaTx = np.array(thetaTx)
        else:
            thetaTx = np.array([thetaTx])
    if not isinstance(thetaRx, np.ndarray):
        if (numThetaRx > 1):
            thetaRx = np.array(thetaRx)
        else:
            thetaRx = np.array([thetaRx])
    
    thetaEff = sqrt(2) * thetaTx * thetaRx / (sqrt(thetaTx ** 2 + thetaRx ** 2)) # Eqn 3.20
    
    return thetaEff

def refractiveidx(h, VaporDensity = 7.5, ScaleHeight = 2000, LatitudeModel = "Standard", Season = "Summer", AtmosphereMeasurements = None):
    '''
    % refractiveidx    Calculates the refractive index
    %   RIDX = refractiveidx(H) calculates the refractive index RIDX using the
    %   International Telecommunication Union (ITU) standard atmospheric model
    %   known as the Mean Annual Global Reference Atmosphere (MAGRA), which
    %   approximates the U.S. Standard Atmosphere 1976 with insignificant
    %   relative error. The input H is an M-length vector of the geometric
    %   heights (altitude above mean sea level [MSL]) in m. The output RIDX is
    %   the refractive index as an M-length row-vector.
    %
    %   The ITU model is valid for altitudes in the range of 0 - 100 km.
    %   Outside of this validity region, the function outputs NaNs.
    %
    %   RIDX = refractiveidx(...,'WaterVaporDensity',RHO0) specifies the
    %   standard ground-level water vapor density RHO0 as a nonnegative scalar
    %   in g/m^3. Applicable only for the default standard model (Mean Annual
    %   Global Reference Atmosphere). Defaults to 7.5 g/m^3.
    %
    %   RIDX = refractiveidx(...,'ScaleHeight',H0) specifies the scale height
    %   H0 (altitude above MSL) as a nonnegative scalar in m. Applicable only
    %   for the default standard model (Mean Annual Global Reference
    %   Atmosphere). Defaults to 2e3 m. For a dry atmosphere, set H0 to 6e3 m.
    %
    %   RIDX = refractiveidx(...,'LatitudeModel',MODEL) specifies the reference
    %   latitude model used. Defaults to 'Standard'. Applicable string inputs
    %   are as follows:
    %      - 'Standard' 
    %           This model is the Mean Annual Global Reference Atmosphere
    %           (MAGRA) that reflects the mean annual temperature and pressure
    %           averaged across the world. Default model.
    %      - 'Low' 
    %           This model is for low latitudes less than 22 deg, where there
    %           exists little seasonal variation.
    %      - 'Mid' 
    %           This model is for mid-latitudes between 22 deg and 45 deg with
    %           seasonal profiles for 'Summer' and 'Winter', which can be
    %           specified using the 'Season' name-value pair as discussed more
    %           below. Defaults to 'Summer' profile.
    %      - 'High' 
    %           This model is for high-latitudes greater than 45 deg with
    %           seasonal profiles for 'Summer' and 'Winter', which can be
    %           specified using the 'Season' name-value pair as discussed more
    %           below. Defaults to 'Summer' profile.
    %   The ITU models are valid for altitudes in the range of 0 - 100 km.
    %   Outside of this validity region, the function outputs NaNs.
    %
    %   RIDX = refractiveidx(...,'Season',SEASON) specifies the season as a
    %   string. Valid only for the 'Mid' and 'High' latitude models  other
    %   models will ignore this input. Applicable inputs are 'Summer' or
    %   'Winter'. Defaults to 'Summer'.
    %
    %   RIDX = refractiveidx(...,'AtmosphereMeasurements',TPWVDENH) specifies
    %   custom atmospheric measurements for the calculation of the refractive
    %   index. TPWVDENH is an N x 4 matrix, where N corresponds to the number
    %   of altitude measurements. The first column is the atmospheric
    %   temperature in K, the second column is the atmospheric dry air pressure
    %   in hPa, the third column is the water vapor density in g/m^3, and the
    %   fourth column is the MSL altitude of the measurements in m. When a
    %   custom model is provided, all other name-value pair options are ignored
    %   and the output refractive index is applicable for the input height H.
    %   If the requested heights H are outside of the range of altitudes
    %   specified within the measurements TPWVDENH, the function outputs NaNs.
    %
    %   [RIDX,N] = refractiveidx(...) additionally outputs the refractivity N
    %   as an M-length row-vector.
    %
    %   refractiveidx(...) with no output arguments, plots the refractive index
    %   n as a function of altitude in km.
    %
    %   % Example 1:
    %   %   Plot the refractive index profile n using the standard global
    %   %   reference model.
    %   h = 0:10:100e3  % m
    %   refractiveidx(h)
    %  
    %   % Example 2: 
    %   %   Compute the refractive index n and the refractivity N at a height 
    %   %   of 20 km using the Mid-Latitude model during winter.
    %   h = 20e3  % m
    %   [n,N] = refractiveidx(h,'LatitudeModel','Mid','Season','Winter')
    %
    %   See also atmositu, tropopl, lenspl, effearthradius.
    
    %   Copyright 2020-2023 The MathWorks, Inc.
    
    %   Reference
    %   [1] International Telecommunication Union (ITU). "The Radio Refractive
    %       Index: Its Formula and Refractivity Data." Recommendation ITU-R
    %       P.453-11, P Series, Radiowave Propagation, July 2015.
    '''

    if isscalar(h):
        h = np.array([h])
    elif not isinstance(h, np.ndarray):
        h = np.array(h)
        
    hgeomKm = h * 1e-3 # Geometric height (km)
    T, P, wvden = stdatm(hgeomKm, VaporDensity, ScaleHeight)

    if AtmosphereMeasurements is None:
        tpwvdenh = [singcolvec(T), singcolvec(P), singcolvec(wvden), singcolvec(h)]
    else:
        tpwvdenh = AtmosphereMeasurements
    modelStr = getModelStr(LatitudeModel, Season)
        
    # Obtain temperature, pressure, and density
    numh = numel(h)
    T = Nan(numh)
    P = Nan(numh)
    rho = Nan(numh)
    
    match modelStr.lower():
        case 'custom':
            idx = np.where((h >= tpwvdenh[0, 3]) & (h <= tpwvdenh[-1, 3]))
            T[idx] = interp1(tpwvdenh[:, 3], tpwvdenh[:, 0], h[idx])
            P[idx] = interp1(tpwvdenh[:, 3],tpwvdenh[:, 1], h[idx])
            rho[idx] = interp1(tpwvdenh[:, 3], tpwvdenh[:, 2], h[idx])
        case _:
            T, P, rho = atmositu(h, VaporDensity = VaporDensity, ScaleHeight = ScaleHeight, LatitudeModel = LatitudeModel, Season = Season)

    ridx, N = ptwvden2RefractiveIndex(T, P, rho) 
    
    return ridx, N

def stdatm(hgeom, rho0, h0):
    '''
    %   This function is for internal use only. It may be removed in the
    %   future.
    
    % stdatm   Standard atmosphere model
    %   [T,P,RHO] = stdatm(HGEOM,RHO0,H0) calculates the ITU standard
    %   atmospheric model based on global averages for temperature, pressure,
    %   and water vapor density. The input HGEOM is an M-length row-vector of
    %   the geometric heights in km, RHO0 is the standard ground-level
    %   water-vapor density in g/m^3 (default is 7.5 g/m^3), and H0 is the
    %   scale height (altitude above MSL) in km (default is 2 km). The output T
    %   is the temperature in K, P is the pressure in hPa, and RHO is the
    %   water-vapor density in g/m^3. All outputs are 1xM length row vectors. 
    
    %   Copyright 2020-2022 The MathWorks, Inc.
    
    %   Reference 
    %   [1] International Telecommunication Union (ITU). "Reference Standard
    %       Atmospheres." Recommendation ITU-R P.835-6, P Series, Radiowave
    %       Propagation, Dec. 2017.
    '''
    
    hgeop = 6356.766 * hgeom / (6356.766 + hgeom) # Geopotential height (km)

    if isscalar(hgeop):
        hgeop = np.array([hgeop])

    T = Nan(numel(hgeop))
    P = Nan(numel(hgeop))
    
    # Standard Temperature (K) and Pressure (hPa), First Height Regime
    # h <= 11
    idx = np.where((hgeop <= 11) & (hgeop >= 0))
    if np.size(idx) != 0:
        T[idx] = 288.15 - 6.5 * hgeop[idx]
        P[idx] = 1013.25 * (288.15 / (288.15-6.5 * hgeop[idx])) ** (-34.1632 / 6.5)
    
    # 11 < h <= 20
    idx = np.where((hgeop <= 20) & (hgeop > 11))
    if np.size(idx) != 0:
        T[idx] = 216.65
        P[idx] = 226.3226 * exp(-34.1632 * (hgeop[idx] - 11) / 216.65)
    
    # 20 < h <= 32
    idx = np.where((hgeop <= 32) & (hgeop > 20))
    if np.size(idx) != 0:
        T[idx] = 216.65 + (hgeop[idx] - 20)
        P[idx] = 54.74980 * (216.65 / (216.65 + (hgeop[idx] - 20))) ** 34.1632
    
    # 32 < h <= 47
    idx = np.where((hgeop <= 47) & (hgeop > 32))
    if np.size(idx) != 0:
        T[idx] = 228.65 + 2.8 * (hgeop[idx] - 32)
        P[idx] = 8.680422 * (228.65 / (228.65 + 2.8 * (hgeop[idx] - 32))) ** (34.1632 / 2.8)
    
    # 47 < h <= 51
    idx = np.where((hgeop <= 51) & (hgeop > 47))
    if np.size(idx) != 0:
        T[idx] = 270.65
        P[idx] = 1.109106 * exp(-34.1632 * (hgeop[idx] - 47) / 270.65)
    
    # 51 < h <= 71
    idx = np.where((hgeop <= 71) & (hgeop > 51))
    if np.size(idx) != 0:
        T[idx] = 270.65 - 2.8 * (hgeop[idx] - 51)
        P[idx] = 0.6694167 * (270.65 / (270.65 - 2.8 * (hgeop[idx] - 51))) ** (-34.1632 / 2.8)
    
    # 71 < h <= 84.852
    idx = np.where((hgeop <= 84.852) & (hgeop > 71))
    if np.size(idx) != 0:
        T[idx] = 214.65 - 2.0 * (hgeop[idx] - 71)
        P[idx] = 0.03956649 * (214.65 / (214.65 - 2.0 * (hgeop[idx] - 71))) ** (-34.1632 / 2.0)
    
    # Standard Temperature (K), Second Height Regime
    # 86 =< h <= 91
    idx = np.where((hgeom <= 91) & (hgeom >= 86))
    if np.size(idx) != 0:
        T[idx] = 186.8673
    
    # 91 < h <= 100
    idx = np.where((hgeom <= 100) & (hgeom > 91))
    if np.size(idx) != 0:
        T[idx] = 263.1905 - 76.3232 * sqrt(1 - ((hgeom[idx] - 91) / 19.9429) ** 2)
    
    # Standard Pressure (hPa), Second Height Regime
    idx = np.where((hgeom <= 100) & (hgeom >= 86))
    if np.size(idx) != 0:
        a0 = 95.571899
        a1 = -4.011801
        a2 = 6.424731e-2
        a3 = -4.789660e-4
        a4 = 1.340543e-6
        P[idx] = exp(a0 + a1 * hgeom[idx] + a2 * hgeom[idx] ** 2 + a3 * hgeom[idx] ** 3 + a4 * hgeom[idx] ** 4)
    
    # Water-vapor density (g/m^3)
    h0 = h0 * 1e-3
    rho = rho0 * exp(-hgeom / h0)
    
    # Calculate water vapor pressure
    e = rho * T / 216.7 # hPa
    
    # Water-vapor density decreases exponentially with increasing altitude up
    # to an altitude where the mixing ratio e/P = 2e-6. Above this altitude,
    # the mixing ratio is constant.
    if rho0 > sqrt(eps()): # Check for dry atmosphere edge case
        mixingRatioConst = e / P < 2e-6
        if any(mixingRatioConst):
            e[mixingRatioConst] = 2e-6 * P[mixingRatioConst]
        rho = e * 216.7 / T
    
    # Edge case where scale height = 0
    if h0 <= sqrt(eps()):
        rho = zeros(np.size(hgeom))

    if len(T) == 1:
        T = T[0]
    if len(P) == 1:
        P = P[0]
    if len(rho) == 1:
        rho = rho[0]

    return T, P, rho

def getModelStr(model,season):
    '''
    %   This function is for internal use only. It may be removed in the
    %   future.
    
    % getModelStr   Determines model name used for plots
    
    %   Copyright 2020-2023 The MathWorks, Inc.
    '''
    # Define string for plots 
    match model.lower():
        case ('standard', 'custom', 'low'):
            modelStr = model
        case _:
            modelStr = f"{season} {model}"
            
    return modelStr

def atmositu(h, VaporDensity = 7.5, ScaleHeight = 2000, LatitudeModel = "Standard", Season = "Summer"):
    '''
    % atmositu    ITU reference atmospheres 
    %   [T,P,WVDEN] = atmositu(H) calculates the International
    %   Telecommunication Union (ITU) standard atmospheric model known as the
    %   Mean Annual Global Reference Atmosphere (MAGRA), which approximates the
    %   U.S. Standard Atmosphere 1976 with insignificant relative error. The
    %   input H is an M-length vector of the geometric heights (altitude above
    %   mean sea level [MSL]) in m. The output T is the temperature in K, P is
    %   the dry air pressure in hPa, and WVDEN is the water-vapor density in
    %   (g/m^3). The outputs are all M-length row-vectors.
    %
    %   The ITU model is valid for altitudes in the range of 0 - 100 km.
    %   Outside of this validity region, the function outputs NaNs.
    %   
    %   [T,P,WVDEN] = atmositu(...,'WaterVaporDensity',RHO0) specifies the
    %   standard ground-level water vapor density RHO0 as a nonnegative scalar
    %   in g/m^3. Applicable only for the water density calculation in the
    %   default standard model (Mean Annual Global Reference Atmosphere,
    %   Equation 6 in ITU-R P.835-6). Defaults to 7.5 g/m^3.
    %
    %   [T,P,WVDEN] = atmositu(...,'ScaleHeight',H0) specifies the scale height
    %   H0 (altitude above MSL) as a nonnegative scalar in m. Applicable only
    %   for the water density calculation in the default standard model (Mean
    %   Annual Global Reference Atmosphere, Equation 6 in ITU-R P.835-6).
    %   Defaults to 2e3 m. For a dry atmosphere, set H0 to 6e3 m.
    %
    %   [T,P,WVDEN] = atmositu(...,'LatitudeModel',MODEL) specifies the
    %   reference latitude model used. Defaults to 'Standard'. Applicable
    %   string inputs are as follows:
    %      - 'Standard' 
    %           This model is the Mean Annual Global Reference Atmosphere
    %           (MAGRA) that reflects the mean annual temperature and pressure
    %           averaged across the world. Default model.
    %      - 'Low' 
    %           This model is for low latitudes less than 22 deg, where there
    %           exists little seasonal variation.
    %      - 'Mid' 
    %           This model is for mid-latitudes between 22 deg and 45 deg with
    %           seasonal profiles for 'Summer' and 'Winter', which can be
    %           specified using the 'Season' name-value pair as discussed more
    %           below. Defaults to 'Summer' profile.
    %      - 'High' 
    %           This model is for high-latitudes greater than 45 deg with
    %           seasonal profiles for 'Summer' and 'Winter', which can be
    %           specified using the 'Season' name-value pair as discussed more
    %           below. Defaults to 'Summer' profile.
    %   The ITU models are valid for altitudes in the range of 0 - 100 km.
    %   Outside of this validity region, the function outputs NaNs.
    %
    %   [T,P,WVDEN] = atmositu(...,'Season',SEASON) specifies the season as a
    %   string. Valid only for the 'Mid' and 'High' latitude models  other
    %   models will ignore this input. Applicable inputs are 'Summer' or
    %   'Winter'. Defaults to 'Summer'.
    %
    %   atmositu(...) with no output arguments, plots the following:
    %      - Temperature (K) versus altitude (km), linear scale
    %      - Pressure (hPa) versus altitude (km), logarithmic x-scale
    %      - Water vapor density (g/m^3) versus altitude (km), logarithmic
    %        x-scale
    %
    %   % Example 1:
    %   %   Plot the temperature, pressure, and water vapor density profiles 
    %   %   using the standard global reference model. 
    %   h = 0:10:100e3  % m
    %   atmositu(h)
    %  
    %   % Example 2: 
    %   %   Compute the temperature, pressure, and water vapor density at a 
    %   %   height of 30 km using the Mid-Latitude model during winter. 
    %   h = 30e3  % m
    %   [T,P,wvden] = atmositu(h,'LatitudeModel','Mid','Season','Winter')
    %
    %   See also refractiveidx, tropopl, lenspl.
    
    %   Copyright 2020-2023 The MathWorks, Inc.
    
    %   Reference
    %   [1] International Telecommunication Union (ITU). "Reference Standard
    %       Atmospheres." Recommendation ITU-R P.835-6, P Series, Radiowave
    %       Propagation, Dec. 2017.
    '''
        
    # Define geometric and geopotential height
    hgeomKm = h * 1e-3 # Geometric height (km)
    
    # Calculate temperature (K), pressure (hPa), and water-vapor density
    # (g/m^3)
    numh = numel(h) 
    T = Nan(numh) 
    P = Nan(numh)
    wvden = Nan(numh)
    
    match LatitudeModel.lower():
        case 'low':
            T, P, wvden = lowatm(hgeomKm) 
        case 'mid':
            match Season.lower():
                case 'summer':
                    T, P, wvden = summermidatm(hgeomKm)
                case 'winter':
                    T, P, wvden = wintermidatm(hgeomKm)
        case 'high':
            match Season.lower():
                case 'summer':
                    T, P, wvden = summerhighatm(hgeomKm)
                case 'winter':
                    T, P, wvden = winterhighatm(hgeomKm)
        case _: # 'standard'
            T, P, wvden = stdatm(hgeomKm, VaporDensity, ScaleHeight)

    return T, P, wvden

def lowatm(h):
    '''
    %   This function is for internal use only. It may be removed in the
    %   future.
    
    % lowatm    Low-latitude annual reference atmosphere model
    %   [T,P,RHO] = lowatm(H) calculates the ITU model for pressure,
    %   temperature, and water vapor density for low-latitudes less that 22
    %   degrees with little seasonal variation. The input H is an M-length
    %   row-vector of the geometric heights in km. The output T is the
    %   temperature in K, P is the pressure in hPa, and RHO is the water-vapor
    %   density in g/m^3. All outputs are 1xM length row vectors. 
    
    %   Copyright 2020 The MathWorks, Inc.
    
    %   Reference 
    %   [1] International Telecommunication Union (ITU). "Reference Standard
    %       Atmospheres." Recommendation ITU-R P.835-6, P Series, Radiowave
    %       Propagation, Dec. 2017.
    '''
    
    T = Nan(numel(h))
    P = Nan(numel(h))
    rho = Nan(numel(h))
    
    # Temperature (K)
    idx = np.where((h < 17) & (h >= 0))
    if np.size(idx) != 0:
        T[idx] = 300.4222 - 6.3533 * h[idx] + 0.005886 * h[idx] ** 2
    
    idx = np.where((h < 47) & (h >= 17))
    if np.size(idx) != 0:
        T[idx] = 194 + (h[idx] - 17) * 2.533
    
    idx = np.where((h < 52) & (h >= 47))
    if np.size(idx) != 0:
        T[idx] = 270
    
    idx = np.where((h < 80) & (h >= 52))
    if np.size(idx) != 0:
        T[idx] = 270 - (h[idx] - 52) * 3.0714
    
    idx = np.where((h <= 100) & (h >= 80))
    if np.size(idx) != 0:
        T[idx] = 184
    
    # Pressure (hPa)
    idx = np.where((h <= 10) & (h >= 0))
    if np.size(idx) != 0:
        P[idx] = 1012.0306 - 109.0338 * h[idx] + 3.6316 * h[idx] ** 2
    
    idx = np.where((h <= 72) & (h > 10))
    if np.size(idx) != 0:
        P10 = 1012.0306 - 109.0338 * 10 + 3.6316 * 10 ** 2 
        P[idx] = P10 * exp(-0.147 * (h[idx] - 10))
    
    idx = np.where((h <= 100) & (h > 72))
    if np.size(idx) != 0:
        P72 = P10 * exp(-0.147 * (72 - 10))
        P[idx] = P72 * exp(-0.165 * (h[idx] - 72))
    
    # Water vapor density (g/m^3)
    idx = np.where((h <= 15) & (h >= 0))
    if np.size(idx) != 0:
        rho[idx] = 19.6542 * exp(-0.2313 * h[idx] - 0.1122 * h[idx] ** 2 + 0.01351 * h[idx] ** 3 - 0.0005923 * h[idx] ** 4)
    
    idx = np.where((h <= 100) & (h > 15))
    if np.size(idx) != 0:
        rho[idx] = 0

    if len(T) == 1:
        T = T[0]
    if len(P) == 1:
        P = P[0]
    if len(rho) == 1:
        rho = rho[0]

    return T, P, rho

def summermidatm(h):
    '''
    %   This function is for internal use only. It may be removed in the
    %   future.
    
    % summermidatm   Mid-latitude reference atmosphere model for summer
    %   [T,P,RHO] = summermidatm(H) calculates the ITU summer season model for
    %   pressure, temperature, and water vapor density for mid-latitudes
    %   between 22 and 45 degrees. The input H is an M-length row-vector of the
    %   geometric heights in km. The output T is the temperature in K, P is the
    %   pressure in hPa, and RHO is the water-vapor density in g/m^3. All
    %   outputs are 1xM length row vectors.
    
    %   Copyright 2020 The MathWorks, Inc.
    
    %   Reference 
    %   [1] International Telecommunication Union (ITU). "Reference Standard
    %       Atmospheres." Recommendation ITU-R P.835-6, P Series, Radiowave
    %       Propagation, Dec. 2017.
    '''
    
    T = Nan(numel(h))
    P = Nan(numel(h))
    rho = Nan(numel(h))
    
    # Temperature
    idx = np.where((h < 13) & (h >= 0))
    if np.size(idx) != 0:
        T[idx] = 294.9838 - 5.2159 * h[idx] - 0.07109 * h[idx] ** 2
    
    idx = np.where((h >= 13) & (h < 17))
    if np.size(idx) != 0:
        T[idx] = 215.15
    
    idx = np.where((h >= 17) & (h < 47))
    if np.size(idx) != 0:
        T[idx] = 215.15 * exp((h[idx] - 17) * 0.008128)
    
    idx = np.where((h >= 47) & (h < 53))
    if np.size(idx) != 0:
        T[idx] = 275
    
    idx = np.where((h >= 53) & (h < 80))
    if np.size(idx) != 0:
        T[idx] = 275 + (1 - exp((h[idx] - 53) * 0.06)) * 20
    
    idx = np.where((h >= 80) & (h <= 100))
    if np.size(idx) != 0:
        T[idx] = 175
    
    # Pressure
    idx = np.where((h <= 10) & (h >= 0))
    if np.size(idx) != 0:
        P[idx] = 1012.8186 - 111.5569 * h[idx] + 3.8646 * h[idx] ** 2 
    
    idx = np.where((h > 10) & (h <= 72))
    if np.size(idx) != 0:
        P10 = 1012.8186 - 111.5569 * 10 + 3.8646 * 10 ** 2
        P[idx] = P10 * exp(-0.147 * (h[idx] - 10))
    
    idx = np.where((h > 72) & (h <= 100))
    if np.size(idx) != 0:
        P72 = P10 * exp(-0.147 * (72 - 10))  
        P[idx] = P72 * exp(-0.165 * (h[idx] - 72))  
    
    # Water vapor density 
    idx = np.where((h <= 15) & (h >= 0))
    if np.size(idx) != 0:
        rho[idx] = 14.3542 * exp(-0.4174 * h[idx] - 0.02290 * h[idx] ** 2 + 0.001007 * h[idx] ** 3)  
    
    idx = np.where((h <= 100) & (h > 15))
    if np.size(idx) != 0:
        rho[idx] = 0 

    if len(T) == 1:
        T = T[0]
    if len(P) == 1:
        P = P[0]
    if len(rho) == 1:
        rho = rho[0]
        
    return T, P, rho

def wintermidatm(h):
    '''
    %   This function is for internal use only. It may be removed in the
    %   future.
    
    % wintermidatm   Mid-latitude reference atmosphere model for winter
    %   [T,P,RHO] = wintermidatm(H) calculates the ITU winter season model for
    %   pressure, temperature, and water vapor density for mid-latitudes
    %   between 22 and 45 degrees. The input H is an M-length row-vector of the
    %   geometric heights in km. The output T is the temperature in K, P is the
    %   pressure in hPa, and RHO is the water-vapor density in g/m^3. All
    %   outputs are 1xM length row vectors.
    
    %   Copyright 2020 The MathWorks, Inc.
    
    %   Reference 
    %   [1] International Telecommunication Union (ITU). "Reference Standard
    %       Atmospheres." Recommendation ITU-R P.835-6, P Series, Radiowave
    %       Propagation, Dec. 2017.
    '''
    
    T = Nan(numel(h))  
    P = Nan(numel(h))  
    rho = Nan(numel(h))
    
    # Temperature
    idx = np.where((h < 10) & (h >= 0))
    if np.size(idx) != 0:
        T[idx] = 272.7241 - 3.6217 * h[idx] - 0.1759 * h[idx] ** 2 
    
    idx = np.where((h >= 10) & (h < 33))
    if np.size(idx) != 0:
        T[idx] = 218 
    
    idx = np.where((h >= 33) & (h < 47))
    if np.size(idx) != 0:
        T[idx] = 218 + (h[idx] - 33) * 3.3571 
    
    idx = np.where((h >= 47) & (h < 53))
    if np.size(idx) != 0:
        T[idx] = 265 
    
    idx = np.where((h >= 53) & (h < 80))
    if np.size(idx) != 0:
        T[idx] = 265 - (h[idx] - 53) * 2.0370 
    
    idx = np.where((h >= 80) & (h <= 100))
    if np.size(idx) != 0:
        T[idx] = 210 
    
    # Pressure
    idx = np.where((h <= 10) & (h >= 0))
    if np.size(idx) != 0:
        P[idx] = 1018.8627 - 124.2954 * h[idx] + 4.8307 * h[idx] ** 2  
    
    idx = np.where((h > 10) & (h <= 72))
    if np.size(idx) != 0:
        P10 = 1018.8627 - 124.2954 * 10 + 4.8307 * 10 ** 2  
        P[idx] = P10 * exp(-0.147 * (h[idx] - 10)) 
    
    idx = np.where((h > 72) & (h <= 100))
    if np.size(idx) != 0:
        P72 = P10 * exp(-0.147 * (72 - 10))
        P[idx] = P72 * exp(-0.155 * (h[idx] - 72))  
    
    # Water vapor density 
    idx = np.where((h <= 10) & (h >= 0))
    if np.size(idx) != 0:
        rho[idx] = 3.4742 * exp(-0.2697 * h[idx] - 0.03604 * h[idx] ** 2 + 0.0004489 * h[idx] ** 3) 
    
    idx = np.where((h <= 100) & (h > 10))
    if np.size(idx) != 0:
        rho[idx] = 0 

    if len(T) == 1:
        T = T[0]
    if len(P) == 1:
        P = P[0]
    if len(rho) == 1:
        rho = rho[0]
        
    return T, P, rho

def summerhighatm(h):
    '''
    %   This function is for internal use only. It may be removed in the
    %   future.
    
    % summerhighatm   High-latitude reference atmosphere model for summer
    %   [T,P,RHO] = summerhighatm(H) calculates the ITU summer season model for
    %   pressure, temperature, and water vapor density for high-latitudes
    %   greater than 45 degrees. The input H is an M-length row-vector of the
    %   geometric heights in km. The output T is the temperature in K, P is the
    %   pressure in hPa, and RHO is the water-vapor density in g/m^3. All
    %   outputs are 1xM length row vectors.
    
    %   Copyright 2020 The MathWorks, Inc.
    
    %   Reference 
    %   [1] International Telecommunication Union (ITU). "Reference Standard
    %       Atmospheres." Recommendation ITU-R P.835-6, P Series, Radiowave
    %       Propagation, Dec. 2017.
    '''
    
    T = Nan(numel(h))  
    P = Nan(numel(h))  
    rho = Nan(numel(h))  
    
    # Temperature
    idx = np.where((h < 10) & (h >= 0))
    if np.size(idx) != 0:
        T[idx] = 286.8374 - 4.7805 * h[idx] - 0.1402 * h[idx] ** 2  
    
    idx = np.where((h >= 10) & (h < 23))
    if np.size(idx) != 0:
        T[idx] = 225 
    
    idx = np.where((h >= 23) & (h < 48))
    if np.size(idx) != 0:
        T[idx] = 225 * exp((h[idx] - 23) * 0.008317) 
    
    idx = np.where((h >= 48) & (h < 53))
    if np.size(idx) != 0:
        T[idx] = 277 
    
    idx = np.where((h >= 53) & (h < 79))
    if np.size(idx) != 0:
        T[idx] = 277 - (h[idx] - 53) * 4.0769 
    
    idx = np.where((h >= 79) & (h <= 100))
    if np.size(idx) != 0:
        T[idx] = 171 
    
    # Pressure
    idx = np.where((h <= 10) & (h >= 0))
    if np.size(idx) != 0:
        P[idx] = 1008.0278 - 113.2494 * h[idx] + 3.9408 * h[idx] ** 2 
    
    idx = np.where((h > 10) & (h <= 72))
    if np.size(idx) != 0:
        P10 = 1008.0278 - 113.2494 * 10 + 3.9408 * 10 ** 2 
        P[idx] = P10 * exp(-0.140 * (h[idx] - 10))
    
    idx = np.where((h > 72) & (h <= 100))
    if np.size(idx) != 0:
        P72 = P10 * exp(-0.140 * (72 - 10))  
        P[idx] = P72 * exp(-0.165 * (h[idx] - 72))  
    
    # Water vapor density 
    idx = np.where((h <= 15) & (h >= 0))
    if np.size(idx) != 0:
        rho[idx] = 8.988 * exp(-0.3614 * h[idx] - 0.005402 * h[idx] ** 2 - 0.001955 * h[idx] ** 3) 
    
    idx = np.where((h <= 100) & (h > 15))
    if np.size(idx) != 0:
        rho[idx] = 0 

    if len(T) == 1:
        T = T[0]
    if len(P) == 1:
        P = P[0]
    if len(rho) == 1:
        rho = rho[0]
    
    return T, P, rho

def winterhighatm(h):
    '''
    %   This function is for internal use only. It may be removed in the
    %   future.
    
    % winterhighatm   High-latitude reference atmosphere model for winter
    %   [T,P,RHO] = winterhighatm(H) calculates the ITU winter season model for
    %   pressure, temperature, and water vapor density for high-latitudes
    %   greater than 45 degrees. The input H is an M-length row-vector of the
    %   geometric heights in km. The output T is the temperature in K, P is the
    %   pressure in hPa, and RHO is the water-vapor density in g/m^3. All
    %   outputs are 1xM length row vectors.
    
    %   Copyright 2020 The MathWorks, Inc.
    
    %   Reference 
    %   [1] International Telecommunication Union (ITU). "Reference Standard
    %       Atmospheres." Recommendation ITU-R P.835-6, P Series, Radiowave
    %       Propagation, Dec. 2017.
    '''
    
    T = Nan(numel(h))  
    P = Nan(numel(h))  
    rho = Nan(numel(h))  
    
    # Temperature
    idx = np.where((h < 8.5) & (h >= 0))
    if np.size(idx) != 0:
        T[idx] = 257.4345 + 2.3474 * h[idx] - 1.5479 * h[idx] ** 2 + 0.08473 * h[idx] ** 3  
    
    idx = np.where((h >= 8.5) & (h < 30))
    if np.size(idx) != 0:
        T[idx] = 217.5 
    
    idx = np.where((h >= 30) & (h < 50))
    if np.size(idx) != 0:
        T[idx] = 217.5 + (h[idx] - 30) * 2.125 
    
    idx = np.where((h >= 50) & (h < 54))
    if np.size(idx) != 0:
        T[idx] = 260 
    
    idx = np.where((h >= 54) & (h <= 100))
    if np.size(idx) != 0:
        T[idx] = 260 - (h[idx] - 54) * 1.667 
    
    # Pressure
    idx = np.where((h <= 10) & (h >= 0))
    if np.size(idx) != 0:
        P[idx] = 1010.8828 - 122.2411 * h[idx] + 4.554 * h[idx] ** 2  
    
    idx = np.where((h > 10) & (h <= 72))
    if np.size(idx) != 0:
        P10 = 1010.8828 - 122.2411 * 10 + 4.554 * 10 ** 2  
        P[idx] = P10 * exp(-0.147 * (h[idx] - 10)) 
    
    idx = np.where((h > 72) & (h <= 100))
    if np.size(idx) != 0:
        P72 = P10 * exp(-0.147 * (72 - 10)) 
        P[idx] = P72 * exp(-0.150 * (h[idx] - 72))  
    
    # Water vapor density 
    idx = np.where((h <= 10) & (h >= 0))
    if np.size(idx) != 0:
        rho[idx] = 1.2319 * exp(0.07481 * h[idx] - 0.0981 * h[idx] ** 2 + 0.00281 * h[idx] ** 3) 
    
    idx = np.where((h <= 100) & (h > 10))
    if np.size(idx) != 0:
        rho[idx] = 0 

    if len(T) == 1:
        T = T[0]
    if len(P) == 1:
        P = P[0]
    if len(rho) == 1:
        rho = rho[0]
        
    return T, P, rho

def ptwvden2RefractiveIndex(T, P, rho):
    '''
    %   This function is for internal use only. It may be removed in the
    %   future.
    
    % ptwvden2RefractiveIndex    Pressure, Temperature, and water vapor density
    % to refractive index
    %   RIDX = ptwvden2RefractiveIndex (T,P,rho) converts the input m-length
    %   row vectors T temperature in K, P the dry air pressure in hPa, and
    %   WVDEN the water-vapor density in (g/m^3) to a refractive index RIDX.
    %
    %   [RIDX,N] = refractiveidx(...) additionally outputs the refractivity N
    %   as an M-length row-vector.
    %
    %   See also refractiveidx, atmositu, tropopl, lenspl, effearthradius.
    
    %   Copyright 2023 The MathWorks, Inc.
    
    %   Reference
    %   [1] International Telecommunication Union (ITU). "The Radio Refractive
    %       Index: Its Formula and Refractivity Data." Recommendation ITU-R
    %       P.453-11, P Series, Radiowave Propagation, July 2015.
    '''

    # Calculate water vapor pressure 
    eWVP = rho * T / 216.7 # hPa 
    
    # Calculate refractivity 
    N = 77.6 * P / T + 72 * eWVP / T + 3.75 * 10 ** 5 * eWVP / T ** 2 # Eqn 2 
    
    # Calculate refractive index 
    ridx = 1 + N * 1e-6 # Eqn 1 
    return ridx, N

def refractionexp(Ns):
    '''
    % refractionexp   CRPL exponential reference atmosphere refraction exponent     
    %   REXP = refractionexp(Ns) calculates the refraction exponent or decay
    %   constant of the CRPL exponential reference atmosphere model in units of
    %   1/km. Ns is the M-length refractivity in N-units at the surface.
    %   Refractivity values at the Earth's surface are typically in the range
    %   between 200 and 450 N-units. Standard atmosphere refractivity is
    %   typically assumed to be 313 N-units.
    %
    %   Atmospheric refraction evidences itself as a deviation in an
    %   electromagnetic ray from a straight line due to variation in air
    %   density as a function of height. The CRPL exponential reference
    %   atmosphere defines the change in refractivity with height as a simple
    %   exponential decay. The CRPL exponential reference model is defined as
    %        N = Ns*exp(-REXP*h),
    %   where h is the height in kilometers. The relationship between the
    %   refractivity N and the index of refraction n is given by
    %        N = (n - 1)*1e6.
    %   Thus, the CRPL exponential reference atmosphere provides the index of
    %   refraction versus height as
    %        n(h) = 1 + Ns*1e-6*exp(-REXP*h).
    %
    %   % Examples:
    %
    %   % Example 1:
    %   %   Calculate the refraction exponents for surface refractivities equal
    %   %   to 200, 313, and 450 N-units. 
    %   rexp = refractionexp([200,313,450])
    %
    %   % Example 2:
    %   %   Plot the radar vertical coverage pattern assuming a default sinc
    %   %   antenna pattern. The frequency is 100 MHz, the antenna height is 10 
    %   %   meters, and the range is 100 km. Assume the surface is smooth, the
    %   %   antenna is not tilted, and the transmitted polarization is
    %   %   horizontal.
    %   freq = 100e6; % Frequency (Hz)
    %   anht = 10;    % Antenna height (m)
    %   rfs  = 100;   % Range (km)
    %   
    %   % Set effective Earth radius using the high latitude atmosphere model
    %   % in winter
    %   [nidx,N] = refractiveidx([0 1e3], ...
    %       'LatitudeModel','High','Season','Winter');
    %   RGradient = (nidx(2) - nidx(1))/1e3;
    %   Re = effearthradius(RGradient); % m
    %   
    %   % Calculate vertical coverage pattern
    %   [vcpKm,vcpangles] = radarvcd(freq,rfs,anht, ...
    %       'EffectiveEarthRadius',Re);
    %   
    %   % Calculate the refraction exponent
    %   Ns = N(1); % Surface refractivity (N-units)
    %   rexp = refractionexp(Ns)
    %
    %   % Plot Blake chart
    %   blakechart(vcpKm,vcpangles, ...
    %       'SurfaceRefractivity',Ns,'RefractionExponent',rexp);
    %
    %   See also refractiveidx, range2height, height2range, height2grndrange,
    %   blakechart, radarvcd.
    
    %   Copyright 2021-2022 The MathWorks, Inc.
    
    %   Reference:
    %   [1] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
    %       Atmosphere. Washington, DC: US Gov. Print. Off, 1959.
    %   [2] Dutton, E. J., and Thayer, G. D. Techniques for Computing
    %       Refraction of Radio Waves in the Troposphere. NBS Technical Note 
    %       97, U.S. National Bureau of Standards, 1961.
    '''
    
    DelN = -7.32 * exp(0.005577 * Ns) # Eqn 1 in Reference 1
    rexp = log(Ns / (Ns + DelN)) # Eqn 3 in Reference 1

    if not isinstance(rexp, np.ndarray):
        rexp = np.array(rexp)
    
    rexp[isinf(rexp)] = 0 # Free space
    
    return rexp

def earthSurfacePermittivity(
    mtl,
    fc,
    temp = None,
    salinity = None,
    liqfrac = None,
    avf = None,
    dsd = None,
    lwc = None,
    sandpercent = None,
    claypercent = None,
    sg = None,
    vwc = None,
    bulkdensity = None,
    gwc = None,
):
    # Copyright 2019-2024 The MathWorks, Inc.

    mtlLib = pd.DataFrame({
        "Name": ["pure-water", "sea-water", "pure-ice", "wet-ice", "multi-year-ice", "dry-snow", "wet-snow", "soil", "vegetation"],
        "Fmax": [1e12, 1e12, 1e12, 1e12, 1e11, 1e11, 1e11, 1e12, 1e12],
        "Tmin": [-4, -4, -60, nan, -30, -60, -60, -273.15, -20],
        "Tmax": [40, 40, 0, nan, -2, 0, 0, inf, inf],
    })
    

    mtlNames = list(mtlLib["Name"].values)
    mtlLib.set_index("Name", inplace=True)
    mtlParams = np.array([item for row in mtlLib.to_numpy() for item in row]).reshape(-1, 3)
    
    # Use lower-case material name and set dry-ice = pure-ice
    if strcmpi(mtl, "dry-ice"):
        mtlLower = "pure-ice"
    else:
        mtlLower = mtl.lower()
    
    # Validate material name
    mtlMatchingIdx = mtlNames.index(mtlLower)
    if isempty(mtlMatchingIdx):
        raise ValueError("Invalid material")
        
    # mtlMatchingIdx = mtlMatchingIdx[0]
    
    # Validate frequency
    Fmax = mtlParams[mtlMatchingIdx, 0]
    if fc > Fmax:
        raise ValueError("Fc out of range")
    fcGHz = fc / 1e9
    
    # Validate temperature
    T = temp
    Tmin = mtlParams[mtlMatchingIdx, 1]
    Tmax = mtlParams[mtlMatchingIdx, 2]
    if ~isnan(Tmin):    
        if isfinite(Tmax):
            if (temp < Tmin) or (temp > Tmax):
                raise ValueError("Temp out of range")
        else:
            if temp < Tmin:
                raise ValueError("Temp low")
    
    # Calculate
    complexEpsilon = []
    match mtlLower:
        case 'pure-water':
            if temp is None:
                raise ValueError("Invalid arguments for specified material")
                
            # Calculate permittivity
            complexEpsilon, _, _, _ = getPureWaterPermittivity(fcGHz, temp)
        case 'sea-water':
            if all(x is None for x in (temp, salinity)):
                raise ValueError("Invalid arguments for specified material")
                
            # Calculate permittivity. Sea water is equal to pure water when
            # salinity equals zero.
            if salinity == 0:
                complexEpsilon, _, _, _ = getPureWaterPermittivity(fcGHz, temp)
            else:
                complexEpsilon = getSeaWaterPermittivity(fcGHz, temp, salinity)
        case 'pure-ice':
            if temp is None:
                raise ValueError("Invalid arguments for specified material")
                
            # Calculate permittivity
            complexEpsilon = getPureIcePermittivity(fcGHz, T)
        case 'wet-ice': # T = 0
            if all(x is None for x in (liqfrac)):
                raise ValueError("Invalid arguments for specified material")
                
            F_wc = liqfrac
    
            # Calculate permittivity
            complexEpsilon = getWetIcePermittivity(fcGHz, F_wc)
        case 'multi-year-ice':
            if all(x is None for x in (temp, avf)):
                raise ValueError("Invalid arguments for specified material")
                
            v_a = avf

            # Calculate permittivity. Multi-year ice is equal to pure ice when
            # the air volume fraction is zero.
            if v_a == 0:
                complexEpsilon = getPureIcePermittivity(fcGHz, T)
            else:
                complexEpsilon = getMultiYearIcePermittivity(fcGHz, T, v_a)
        case 'dry-snow':
            if all(x is None for x in (temp, dsd)):
                raise ValueError("Invalid arguments for specified material")
                
            rho_ds = dsd

            # Calculate permittivity
            complexEpsilon = getDrySnowPermittivity(fcGHz, T, rho_ds)
        case 'wet-snow':
            if all(x is None for x in (temp, dsd, lwc)):
                raise ValueError("Invalid arguments for specified material")
                
            rho_ds = dsd
            F_wc = lwc

            # Calculate permittivity
            complexEpsilon = getWetSnowPermittivity(fcGHz, T, rho_ds, F_wc)
        case 'soil':
            if all(x is None for x in (temp, sandpercent, claypercent, sg, vwc)):
                raise ValueError("Invalid arguments for specified material")
                
            # Validate sand and clay percentage sum no larger than 1
            if claypercent + sandpercent > 100:
                raise ValueError("Invalid soil percentage")

            rho_s = sg
            m_v = vwc

            # Get bulk density
            if bulkdensity is not None:
                rho_b = bulkdensity
            else: # Calculate bulk density using [ITU-R P.527-6, eq. (57)]
                # Ignore a constituent it its percentage if < 1%, but sum
                # up the included consitutent percentages to 100% [ITU-R
                # P.527-6, p. 18]
                siltPercent = (100 - sandpercent - claypercent)
                P = [[sandpercent], [claypercent], [siltPercent]]
                P[P < 1] = 0
                P = P * (100 / sum(P))
                P[P < 1] = 1 # log(1) = 0

                mult = [[0.078886], [0.038753], [0.032732]]
                rho_b = 1.07256 + sum(mult * log(P))

            # Calculate permittivity
            complexEpsilon = getSoilPermittivity(fcGHz, T, sandpercent, claypercent, rho_s, m_v, rho_b)
        case _: #  'vegetation'
            if all(x is None for x in (temp, gwc)):
                raise ValueError("Invalid arguments for specified material")
                
            # Validate gravimetric water content
            M_g = gwc

            # Calculate permittivity depending on temperature
            if T > 0:
                complexEpsilon = getVegetationPermittivityAboveFreezing(fcGHz, T, M_g)
            else:
                complexEpsilon = getVegetationPermittivityBelowFreezing(fcGHz, T, M_g)
    
    # Real relative permittivity
    epsilon = real(complexEpsilon)
    
    # Absolute permittivity for free space (electric constant)
    epsilon0 = 8.854187817e-12 # From [1]
    
    # Derive conductivity from complex relative permittivity
    sigma = -imag(complexEpsilon) * (2 * pi * fc * epsilon0)
    
    # Conductivity should be non-negative
    if sigma < 0:
        raise ValueError("Invalid conductivity")
    
    return epsilon, sigma, complexEpsilon

def getWaterParams(T):
    # ITU-R P.527-6, Section 5.1.1

    theta = 300 / (T + 273.15) - 1               # Eq. (11)        
    eps_s = 77.66 + 103.3 * theta                # Eq. (8) 
    eps_1 = 0.0671 * eps_s                       # Eq. (9)
    eps_inf = 3.52 - 7.52 * theta                # Eq. (10)
    f_1 = 20.20 - 146.4 * theta + 316 * theta ** 2    # Eq. (12)
    f_2 = 39.8 * f_1                             # Eq. (13)
    
    return eps_s, eps_1, eps_inf, f_1, f_2
    
def getPureWaterPermittivity(fcGHz, T):
    # ITU-R P.527-6, Section 5.1.1
    
    eps_s, eps_1, eps_inf, f_1, f_2 = getWaterParams(T)
    
    # Complex relative permittivity
    eps_r = (eps_s - eps_1)   / (1 + (fcGHz / f_1) ** 2) + (eps_1 - eps_inf) / (1 + (fcGHz / f_2) ** 2) + eps_inf # Eq. (6)
    eps_c = (fcGHz / f_1) * (eps_s - eps_1) / (1 + (fcGHz / f_1) ** 2) + (fcGHz / f_2) * (eps_1 - eps_inf) / (1 + (fcGHz / f_2) ** 2) # Eq. (7)
    eps = eps_r - 1j * eps_c # Eq. (5)
    
    return eps, eps_r, eps_c, f_1
    
def getSeaWaterPermittivity(fcGHz, T, S):
    # ITU-R P.527-6, Section 5.1.2
    
    # Supplemental terms
    eps_s, eps_1, eps_inf, f_1, f_2 = getWaterParams(T)
    eps_ss = eps_s * exp(-3.33330e-3 * S + 4.74868e-6 * S ** 2) # Eq. (17)
    f_1s = f_1 * (1 + S * (2.3232e-3 - 7.9208e-5 * T + 3.6764e-6 * T ** 2 + 3.5594e-7 * T ** 3 + 8.9795e-9 * T ** 4)) # Eq. (18) in GHz
    eps_1s = eps_1 * exp(-6.28908e-3 * S + 1.76032e-4 * S ** 2 - 9.22144e-5 * T * S) # Eq. (19)
    f_2s = f_2 * (1 + S * (-1.99723e-2 + 1.81176e-4 * T)) # Eq. (20) in GHz
    eps_infs = eps_inf * (1 + S * (-2.04265e-3 + 1.57883e-4 * T)) # Eq. (21)
    sigma_35 = 2.903602 + 8.607e-2 * T + 4.738817e-4 * T ** 2 - 2.991e-6 * T ** 3 + 4.3047e-9 * T ** 4 # Eq. (23)
    R_15 = S * (37.5109 + 5.45216 * S + 1.4409e-2 * S ** 2) / (1004.75 + 182.283 * S + S ** 2) # Eq. (24)
    alpha_0 = (6.9431 + 3.2841 * S - 9.9486e-2 * S ** 2) / (84.850 + 69.024 * S + S ** 2) # Eq. (26)
    alpha_1 = 49.843 - 0.2276 * S + 0.198e-2 * S ** 2 # Eq. (27)
    R_T15 = 1 + alpha_0 * (T - 15) / (alpha_1 + T) # Eq. (25)
    sigma_sw = sigma_35 * R_15 * R_T15 # Eq. (22) in S/m
    
    # Complex relative permittivity
    eps_r = (eps_ss - eps_1s) / (1 + (fcGHz / f_1s) ** 2) + (eps_1s - eps_infs) / (1 + (fcGHz / f_2s) ** 2) + eps_infs # Eq. (15)
    eps_c = (fcGHz / f_1s) * (eps_ss - eps_1s) / (1 + (fcGHz / f_1s) ** 2) + (fcGHz / f_2s) * (eps_1s - eps_infs) / (1 + (fcGHz/f_2s) ** 2) + 18 * sigma_sw / fcGHz # Eq. (16)
    eps = eps_r - 1j * eps_c # Eq. (14)
    
    return eps
    
def getPureIcePermittivity(fcGHz, T):
    # ITU-R P.527-6, Section 5.1.3.1
    
    # Supplemental terms
    theta = 300 / (T + 273.15) - 1           # Eq. (34)  
    A = (0.00504 + 0.0062 * theta) * exp(-22.1 * theta) # Eq. (31)
    tau = 335 / (T + 273.15)               # Eq. (33)
    B = 0.0207 / (T + 273.15) * exp(-tau) / (exp(-tau) - 1) ** 2 + 1.16e-11 * fcGHz ** 2 + exp(-9.963 + 0.0372 * T) # Eq. (32)
    
    # Complex relative permittivity
    eps_r = 3.1884 + 0.00091 * T             # Eq. (29)
    eps_c = A / fcGHz + B * fcGHz          # Eq. (30)
    eps = eps_r - 1j * eps_c                 # Eq. (28)
    
    return eps
    
def getWetIcePermittivity(fcGHz, F_wc):
    # ITU-R P.527-5, Section 5.1.3.2; not in ITU-R P.527-6
    
    eps_ice = getPureIcePermittivity(fcGHz, 0)
    eps_pw, _, _, _  = getPureWaterPermittivity(fcGHz, 0)
    
    # Complex relative permittivity
    eps = (((eps_ice + 2 * eps_pw) + 2 * (eps_ice - eps_pw) * (1 - F_wc)) / ((eps_ice + 2 * eps_pw) - (eps_ice-eps_pw) * (1 - F_wc))) * eps_pw # Eq. (35)
     
    return eps
    
def getMultiYearIcePermittivity(fcGHz, T, v_a):
    # ITU-R P.527-6, Section 5.1.3.3.2
    
    # Supplemental terms [p. 15]
    eps_ice = getPureIcePermittivity(fcGHz, T)
    
    # Complex relative permittivity [eq. (49)-(50)]
    A = 2
    B = 1 - 2 * eps_ice - (3 * v_a) * (1 - eps_ice)
    C = -eps_ice
    eps = (-B + sqrt(B ** 2 - 4 * A * C)) / (2 * A) # Reference has wrong sign for sqrt
    
    # Imaginary permittivity can exceed zero due to numerical precision
    eps = real(eps) + 1j * matmin(imag(eps), 0)
    
    return eps
    
def getDrySnowPermittivity(fcGHz, T, rho_ds):
    # ITU-R P.527-6, Section 5.1.4.1
    
    # Supplemental terms [p. 15]
    eps_ice = getPureIcePermittivity(fcGHz, T)
    rho_ice = 0.916 # g/cm^3
    f_ice = rho_ds / rho_ice
    
    # Complex relative permittivity [eq. (51)-(53)]
    if rho_ds <= 0.5:
        eps_r = 1 + 1.9 * rho_ds
    else:
        eps_r = 0.51 + 2.88 * rho_ds

    eps_i = -3 * imag(eps_ice) * f_ice * eps_r ** 2 * (2 * eps_r + 1) / ((real(eps_ice) + 2 * eps_r) * (real(eps_ice) + 2 * eps_r ** 2))
    eps = eps_r - 1j * eps_i
    
    return eps
    
def getWetSnowPermittivity(fcGHz, T, rho_ds, F_wc):
    # ITU-R P.527-6, Section 5.1.4.2
    
    # Supplemental terms [p. 16]
    eps_ds = getDrySnowPermittivity(fcGHz, T, rho_ds)
    eps_pw, _, _, _ = getPureWaterPermittivity(fcGHz, T)
    
    # Complex relative permittivity [eq. (54)-(55)]
    A = 2
    B = eps_pw - 2 * eps_ds - 3 * F_wc * (eps_pw - eps_ds)
    C = -eps_pw * eps_ds
    eps = (-B + sqrt(B ** 2 - 4 * A * C)) / (2 * A)
    
    # Imaginary permittivity can exceed zero due to numerical precision
    eps = real(eps) + 1j * matmin(imag(eps), 0)
    
    return eps
    
def getSoilPermittivity(fcGHz, T, P_sand, P_clay, rho_s, m_v, rho_b):
    # ITU-R P.527-6, Section 5.2
    
    # Complex relative permittivity of free water [eq. (65)-(70)]
    sigma_1 = 0.0467 + 0.2204 * rho_b - 0.004111 * P_sand - 0.006614 * P_clay 
    sigma_2 = -1.645 + 1.939 * rho_b - 0.0225622 * P_sand + 0.01594 * P_clay
    sigma_eff_r = (fcGHz / 1.35) * (sigma_1 - sigma_2) / (1 + (fcGHz / 1.35) ** 2)
    sigma_eff_i = sigma_2 + (sigma_1 - sigma_2) / (1 + (fcGHz / 1.35) ** 2)
    _, eps_pw_r, eps_pw_c, _ = getPureWaterPermittivity(fcGHz, T)
    temp = 18 * (rho_s - rho_b) / (fcGHz * rho_s * m_v)
    eps_fw_r = eps_pw_r + temp * sigma_eff_r
    eps_fw_i = eps_pw_c + temp * sigma_eff_i
    
    # Complex relative permittivity [eq. (58)-(64)]
    alpha = 0.65
    beta_r = 1.2748 - 0.00519 * P_sand - 0.00152 * P_clay
    beta_i = 1.33797 - 0.00603 * P_sand - 0.00166 * P_clay
    eps_sm_r = (1.01 + 0.44 * rho_s) ** 2 - 0.062
    eps_r = (1 + (rho_b / rho_s) * (eps_sm_r ** alpha - 1) + (m_v ** beta_r) * (eps_fw_r ** alpha) - m_v) ** (1/alpha)
    eps_i = ((m_v ** beta_i) * (eps_fw_i ** alpha)) ** (1 / alpha)
    eps = eps_r - 1j * eps_i
    
    return eps
    
def getVegetationPermittivityAboveFreezing(fcGHz, T, M_g):
    # ITU-R P.527-6, Section 5.3.1
    
    # Supplemental terms [eq. (75)-(77)]
    eps_dv = 1.7 - 0.74 * M_g + 6.16 * M_g ** 2
    v_fw = M_g * (0.55 * M_g - 0.076)
    v_bw = 4.64 * M_g ** 2 / (1 + 7.36 * M_g ** 2)
    
    # Complex relative permittivity [eq. (72)-(74)]
    _, eps_pw_r, eps_pw_i, f_1 = getPureWaterPermittivity(fcGHz, T)
    temp1 = sqrt(fcGHz / (0.02 * f_1))
    temp2 = fcGHz / (0.01 * f_1)
    eps_r = eps_dv + v_fw * eps_pw_r + v_bw * (2.9 + 55 * (1 + temp1) / (1 + 2 * temp1 + temp2))
    eps_i = v_fw * (eps_pw_i + 22.86 / fcGHz) + v_bw * (55 * temp1 / (1 + 2 * temp1 + temp2))
    eps = eps_r - 1j * eps_i
    
    return eps
    
def getVegetationPermittivityBelowFreezing(fcGHz, T, M_g):
    # ITU-R P.527-6, Section 5.3.2
    
    # Supplemental terms [eq. (80)-(89)]
    T_f = -6.5 # [p. 24]    
    Delta = T - T_f
    c1 = fcGHz / 1.2582
    c2 = 0.2054
    c2a = c2 * pi / 2
    c3 = 0.4108
    denom = 1 + 2 * c1 ** c2 * cos(c2a) + c1 ** c3
    X1 = (1 + c1 ** c2 * cos(c2a)) / denom
    Y1 = c1 ** c2 * sin(c2a) / denom
    A_ice = 0.001 - 0.012 * M_g + 0.0082 * M_g ** 2
    B_ice = 0.036 - 0.2389 * M_g + 0.1435 * M_g ** 2
    C_ice = -0.0538 + 0.4616 * M_g - 0.3398 * M_g ** 2
    v_ice = A_ice * Delta ** 2 + B_ice * Delta + C_ice
    eps_dv = 6.76 - 10.24 * M_g + 6.19 * M_g ** 2
    v_fw = (-0.106 + 0.6591 * M_g - 0.610 * M_g ** 2) * exp((0.06 + 0.6883 * M_g + 0.0001 * M_g ** 2) * Delta)
    v_bw = (-0.16 + 1.1876 * M_g - 0.387 * M_g ** 2) * exp((0.721 - 1.2733 * M_g + 0.8139 * M_g ** 2) * Delta)
    
    # Complex relative permittivity [eq. (72) & (78)-(79)]
    c4 = fcGHz / 9
    c5 = 82.2 / (1 + c4 ** 2)
    eps_r = eps_dv + v_fw * (4.9 + c5) + v_bw * (8.092 + 14.2067 * X1) + 3.15 * v_ice
    eps_i = v_fw * (c4 * c5 + 11.394 / fcGHz) + 14.2067 * v_bw * Y1
    eps = eps_r - 1j * eps_i
    
    return eps

def landroughness(landType):
    '''
    landroughness Surface height standard deviation for land
       HGTSD = landroughness(LANDTYPE) outputs the standard deviation of
       surface height for the specified land type as a scalar in m (i.e.,
       surface roughness).
    
       The input LANDTYPE is a character array specified as one of the
       following:
          - 'RuggedMountains'
          - 'Mountains'
          - 'Metropolitan' 
          - 'Urban'
          - 'WoodedHills'
          - 'RollingHills'
          - 'Woods'
          - 'Farm'
          - 'Desert'
          - 'Flatland'
          - 'Smooth'
    
       [HGTSD,BETA0,VEGTYPE] = landroughness(...) additionally outputs the
       scalar slope of the land type BETA0 in degrees and a character array
       indicating the vegetation type. Note that BETA0 is 1.4 times the RMS
       surface slope.
    
       The assumed vegetation of the different land types are as follows:
       Land Type        |        Vegetation Type
       -----------------------------------------------------------------------
       RuggedMountains  |        Trees (dense)
       Mountains        |        Trees (dense)
       Metropolitan     |        None
       Urban            |        None
       WoodedHills      |        Trees (dense)
       RollingHills     |        Brush (dense)
       Woods            |        Trees (dense)
       Farm             |        Grass (thin)
       Desert           |        Grass (thin)
       Flatland         |        Grass (thin)
       Smooth           |        None
    
        Example:
          Obtain the surface height standard deviation assuming an urban land
          type.
       hgtsd = landroughness('Urban')
    
       See also searoughness, landreflectivity, seareflectivity,
       clutterSurfaceRCS, radarpropfactor, radarvcd, blakechart.
    
       Copyright 2020 The MathWorks, Inc.
    
       Reference
       [1] Barton, David K. Radar Equations for Modern Radar. 1st edition.
           Norwood, MA: Artech House, 2013.
       [2] Long, Maurice W. Radar Reflectivity of Land and Sea. 3rd Ed.
           Boston: Artech House, 2001.
       [3] Nathanson, Fred E., J. Patrick Reilly, and Marvin N. Cohen. Radar
           Design Principles. 2nd Ed. Mendham, NJ: SciTech Publishing, 1999.
    '''
    
    # Setup table (From Table 9.2 and Appendix Example 9-1 Surface Clutter)
    match landType:
        case 'Rugged Mountains' | 'RuggedMountains' | 'Mountains':
            hgtsd = 100           # m
            beta0 = rad2deg(0.1)  # deg
            vegType = 'Trees'     # Character array
        case 'Metropolitan':
            hgtsd = 33            # m (assumes average building is ~ 10 stories tall) 
            beta0 = rad2deg(0.1)  # deg
            vegType = 'None'      # Character array
        case 'Urban':
            hgtsd = 10            # m
            beta0 = rad2deg(0.1)  # deg
            vegType = 'None'      # Character array
        case 'Wooded Hills' | 'WoodedHills':
            hgtsd = 10            # m
            beta0 = rad2deg(0.05) # deg
            vegType = 'Trees'     # Character array
        case 'Rolling Hills' | 'RollingHills':
            hgtsd = 10            # m
            beta0 = rad2deg(0.05); # deg
            vegType = 'Brush'     # Character array
        case 'Woods':
            hgtsd = 3             # m
            beta0 = rad2deg(0.03) # deg
            vegType = 'Trees'     # Character array
        case 'Farm' | 'Farmland' | 'Desert':
            hgtsd = 3             # m
            beta0 = rad2deg(0.03) # deg
            vegType = 'Grass'     # Character array
        case 'Flatland':
            hgtsd = 1             # m
            beta0 = rad2deg(0.02) # deg
            vegType = 'Grass'     # Character array
        case _:
            # 'Smooth'
            hgtsd = 0.3           # m
            beta0 = rad2deg(0.01) # deg
            vegType = 'None'      # Character array

    return hgtsd, beta0, vegType

def tropopl(R, f, platformHeight, el, VaporDensity = 7.5, ScaleHeight = 2000, LatitudeModel = "Standard", Season = "Summer", AtmosphereMeasurements = None):
    '''
    % tropopl    Slant-path loss due to atmosphere gaseous absorption
    %   Lgas = tropopl(R,F,H,EL) calculates the path loss L in dB due to
    %   tropospheric refraction using the International Telecommunication Union
    %   (ITU) standard atmospheric model known as the Mean Annual Global
    %   Reference Atmosphere (MAGRA), which approximates the U.S. Standard
    %   Atmosphere 1976 with insignificant relative error.
    % 
    %   R is an M-length vector indicating the slant range in m, F is the
    %   N-length frequency vector in Hz, H is a scalar indicating the MSL
    %   altitude of the radar platform in m, and EL specifies the scalar or
    %   M-length vector of initial elevation angles (in deg) at the radar
    %   platform.
    %
    %   Lgas is an MxN matrix whose columns represent the path loss of each
    %   propagation path under the corresponding frequency.
    %
    %   Lgas = tropopl(...,'WaterVaporDensity',RHO0) specifies the standard
    %   ground-level water vapor density RHO0 as a nonnegative scalar in g/m^3.
    %   Applicable only for the default standard model (Mean Annual Global
    %   Reference Atmosphere). Defaults to 7.5 g/m^3.
    %
    %   Lgas = tropopl(...,'ScaleHeight',H0) specifies the scale height H0
    %   (altitude above MSL) as a positive scalar in m. Applicable only for the
    %   default standard model (Mean Annual Global Reference Atmosphere).
    %   Defaults to 2e3 m. For a dry atmosphere, set H0 to 6e3 m.
    %
    %   Lgas = tropopl(...,'LatitudeModel',MODEL) specifies the reference
    %   latitude model used. Defaults to 'Standard'. Applicable string inputs
    %   are as follows:
    %      - 'Standard' 
    %           This model is the Mean Annual Global Reference Atmosphere
    %           (MAGRA) that reflects the mean annual temperature and pressure
    %           averaged across the world. Default model.
    %      - 'Low' 
    %           This model is for low latitudes less than 22 deg, where there
    %           exists little seasonal variation.
    %      - 'Mid' 
    %           This model is for mid-latitudes between 22 deg and 45 deg with
    %           seasonal profiles for 'Summer' and 'Winter', which can be
    %           specified using the 'Season' name-value pair as discussed more
    %           below. Defaults to 'Summer' profile.
    %      - 'High' 
    %           This model is for high-latitudes greater than 45 deg with
    %           seasonal profiles for 'Summer' and 'Winter', which can be
    %           specified using the 'Season' name-value pair as discussed more
    %           below. Defaults to 'Summer' profile.
    %
    %   Lgas = tropopl(...,'Season',SEASON) specifies the season as a string.
    %   Valid only for the 'Mid' and 'High' latitude models; other models will
    %   ignore this input. Applicable inputs are 'Summer' or 'Winter'. Defaults
    %   to 'Summer'.
    %
    %   Lgas = tropopl(...,'AtmosphereMeasurements',TPWVDENH) specifies custom
    %   atmospheric measurements for the calculation of the attenuation due to
    %   gases. TPWVDENH is an N x 4 matrix, where N corresponds to the number
    %   of altitude measurements. The first column is the atmospheric
    %   temperature in K, the second column is the atmospheric dry air pressure
    %   in hPa, the third column is the water vapor density in g/m^3, and the
    %   fourth and final column is the height of the measurements in m. When a
    %   custom model is provided, all other name-value pair options are
    %   ignored.
    %
    %   The ITU atmosphere gas model is valid between 1 to 1000 GHz.
    %
    %   [Lgas,Llens] = tropopl(...) also outputs the corresponding lens loss L
    %   in dB as an MxN matrix for elevation angles EL less than 50 deg. The
    %   variation in refractivity versus altitude makes the atmosphere act like
    %   a lens with loss independent of frequency. Rays leaving an antenna are
    %   refracted in the troposphere and the energy radiated within some
    %   angular extent is distributed over a slightly greater angular sector,
    %   thereby reducing the energy density relative to propagation in a
    %   vacuum.
    %
    %   Lens loss is not due to a dissipation of energy and does not contribute
    %   to the system noise temperature. Lens loss is a function of range and
    %   elevation angle, and the loss varies slowly at long ranges. Lens loss
    %   can be significant at long ranges and shallow propagation angles.
    %
    %   Examples:
    %
    %   % Example 1: Calculate the zenith attenuation with a standard 
    %   % atmosphere over the range of 1 to 1000 GHz. 
    %   R  = 100e3;             % m
    %   f  = (1:1000).*1e9;     % Hz
    %   ht = 0;                 % m
    %   el = 90;                % deg
    %   Lgas = tropopl(R,f,ht,el);
    % 
    %   % Plot
    %   semilogy(f*1e-9,Lgas);
    %   xlabel('Frequency (GHz)')
    %   ylabel('Attenuation (dB)') 
    %   title('Zenith Attenuation with Standard Atmosphere')
    %   grid on; 
    %   
    %   % Example 2: Calculate the attenuation vs. range for a frequency of 100
    %   % GHz with an elevation of 5 deg using the mid-latitude, winter 
    %   % atmospheric model. 
    %   R  = (10:200)*1e3;      % m
    %   f  = 100e9;             % Hz
    %   ht = 0;                 % m
    %   el = 5;                 % deg
    %   Lgas = tropopl(R,f,ht,el,'LatitudeModel','Mid','Season','Winter');
    % 
    %   % Plot
    %   semilogy(R.*1e-3,Lgas);
    %   xlabel('Range (km)')
    %   ylabel('Attenuation (dB)') 
    %   title('Attenuation for Mid-Latitude, Winter Atmosphere')
    %   grid on; 
    %
    %   See also atmositu, refractiveidx, lenspl, gaspl, effearthradius.
    
    %   Copyright 2020-2022 The MathWorks, Inc.
    
    %   Reference 
    %   [1] International Telecommunication Union (ITU). "Attenuation by
    %       Atmospheric Gases." Recommendation ITU-R P.676-12, P Series,
    %       Radiowave Propagation, Aug. 2019.
    '''

    hgeomKm = platformHeight * 1e-3 # Geometric height (km)
    T, P, wvden = stdatm(hgeomKm, VaporDensity, ScaleHeight)

    if AtmosphereMeasurements is None:
        tpwvdenh = [singcolvec(T), singcolvec(P), singcolvec(wvden), singcolvec(platformHeight)]
    else:
        tpwvdenh = AtmosphereMeasurements
    modelStr = getModelStr(LatitudeModel, Season)
    
    # Earth radius
    Re = Rearth * 1e-3 # km
    
    # Calculate layers 
    platformHeightKm = platformHeight * 1e-3
    edgesITU = calculatelayers(modelStr, tpwvdenh)
    edgesITUm = edgesITU * 1e3
    
    # Check for ducting conditions
    
    # midpointLayers = 0.5 * edgesITUm[1:-1] + 0.5 * edgesITUm[0:-2]
    midpointLayers = 0.5 * edgesITUm[0:-1] + 0.5 * edgesITUm[1:]
    ridx, N = refractiveidx(midpointLayers, VaporDensity = VaporDensity, ScaleHeight = ScaleHeight, LatitudeModel = LatitudeModel, Season = Season, AtmosphereMeasurements = AtmosphereMeasurements)
    refGrad = diff(N) / diff(midpointLayers)
    cond = any(refGrad <= -157) # From ITU-R P.676-11
    if cond:
        raise ValueError("tropopl: refgrad below -157")
    
    # Check for subrefractive conditions
    cond = any(refGrad > 0)
    if cond:
        raise ValueError("tropopl: refgrad cannot be positive")
    
    # Obtain standard pressure, temperature, and water vapor pressure
    match modelStr.lower():
        case 'custom':
            numh = numel(midpointLayers)
            T = Nan(numh)
            P = Nan(numh)
            rho = Nan(numh)
            idx = midpointLayers >= tpwvdenh[0, 3] & midpointLayers <= tpwvdenh[-1, 4]
            if size(tpwvdenh, 0) > 2: # Necessary for codegen. Parser already checks for this. 
                T[idx] = interp1(tpwvdenh[:, 3], tpwvdenh[:, 0], midpointLayers[idx])
                P[idx] = interp1(tpwvdenh[:, 3], tpwvdenh[:, 1], midpointLayers[idx])
                rho[idx] = interp1(tpwvdenh[:, 3], tpwvdenh[:, 2], midpointLayers[idx])
        case _:
            T, P, rho = atmositu(midpointLayers, VaporDensity = VaporDensity, ScaleHeight = ScaleHeight, LatitudeModel = LatitudeModel, Season = Season)

    # Gas and oxygen attenuation loss 
    numF = numel(f)
    gamman = gasatt(f, T, P, rho)
    if np.ndim(gamman) > 1:
        idxNotNaN = ~isnan(gamman[0, :])
    else:
        idxNotNaN = ~isnan(gamman)
    
    numEl = numel(el)
    Lgas = zeros((numel(R), numF))

    # if isscalar(el):
    #     el = [el]

    # for ie in range(0, numEl):
    #     # Convert elevation angle
    #     elRad = el[ie] * pi / 180 # Convert to rad
    
    #     # Propagate
    #     if el[ie] >= 0:
    #         edges, anEdges, alphanEdges, _, hlayers = propagatePosEl(platformHeightKm, elRad, edgesITU, Re, ridx)
    #     else:
    #         edges, anEdges, alphanEdges, _, hlayers = propagateNegEl(platformHeightKm, elRad, edgesITU, Re, ridx)

    #     RsEdges, _, _, _ = refractiongeometry(elRad, anEdges, alphanEdges, edges, Re)
    
    #     Rs = np.hstack((0, RsEdges))
    #     an = np.hstack((0, anEdges))
    #     Rkm, idx = tropopl_getRange(R, ie, numEl)

    #     if numF > 1:
    #         for ifreq in range(0, numF):
    #             # Compute cumulative loss for layers
    #             thisGamman = np.hstack((0, interp1(midpointLayers[idxNotNaN], gamman[ifreq, idxNotNaN], hlayers)))
    #             thisGamman[isnan(thisGamman)] = 0
    #             angamman = cumsum(bsxfun("times", an, thisGamman))
    
    #             # Define Lgas for input ranges
    #             Lgas[idx, ifreq] = interpLosses(Rkm, Rs, angamman)
    #     else:
    #         # Compute cumulative loss for layers
    #         thisGamman = np.hstack((0, interp1(midpointLayers[idxNotNaN], gamman[idxNotNaN], hlayers)))
    #         thisGamman[isnan(thisGamman)] = 0
    #         angamman = cumsum(bsxfun("times", an, thisGamman))

    #         # Define Lgas for input ranges
    #         Lgas = interpLosses(Rkm, Rs, angamman)

    if numEl > 1:
        for ie in range(0, numEl):
            # Convert elevation angle
            elRad = el[ie] * pi / 180 # Convert to rad
        
            # Propagate
            if el[ie] >= 0:
                edges, anEdges, alphanEdges, _, hlayers = propagatePosEl(platformHeightKm,elRad,edgesITU,Re,ridx)
            else:
                edges, anEdges, alphanEdges, _, hlayers = propagateNegEl(platformHeightKm,elRad,edgesITU,Re,ridx);

            RsEdges, _, _, _ = refractiongeometry(elRad, anEdges, alphanEdges, edges, Re)
        
            Rs = np.hstack((0, RsEdges))
            an = np.hstack((0, anEdges))
            Rkm, idx = tropopl_getRange(R, ie, numEl)

            if numF > 1:
                for ifreq in range(0, numF):
                    # Compute cumulative loss for layers
                    thisGamman = np.hstack((0, interp1(midpointLayers[idxNotNaN], gamman[ifreq, idxNotNaN], hlayers)))
                    thisGamman[isnan(thisGamman)] = 0
                    angamman = cumsum(bsxfun("times", an, thisGamman))
            
                    # Define Lgas for input ranges
                    Lgas[idx, ifreq] = interpLosses(Rkm, Rs, angamman)
            else:
                # Compute cumulative loss for layers
                thisGamman = np.hstack((0, interp1(midpointLayers[idxNotNaN], gamman[idxNotNaN], hlayers)))
                thisGamman[isnan(thisGamman)] = 0
                angamman = cumsum(bsxfun("times", an, thisGamman))
        
                # Define Lgas for input ranges
                Lgas[idx, 0] = interpLosses(Rkm, Rs, angamman)
    else:
        # Convert elevation angle
        elRad = el * pi / 180 # Convert to rad
    
        # Propagate
        if el >= 0:
            edges, anEdges, alphanEdges, _, hlayers = propagatePosEl(platformHeightKm, elRad, edgesITU, Re, ridx)
        else:
            edges, anEdges, alphanEdges, _, hlayers = propagateNegEl(platformHeightKm, elRad, edgesITU, Re, ridx)
            
        RsEdges, _, _, _ = refractiongeometry(elRad, anEdges, alphanEdges, edges, Re)

        Rs = np.hstack((0, RsEdges))
        an = np.hstack((0, anEdges))
        Rkm, idx = tropopl_getRange(R, 0, numEl)

        if numF > 1:
            for ifreq in range(0, numF):
                # Compute cumulative loss for layers
                thisGamman = np.hstack((0, interp1(midpointLayers[idxNotNaN], gamman[ifreq, idxNotNaN], hlayers)))
                thisGamman[isnan(thisGamman)] = 0
                angamman = cumsum(bsxfun("times", an, thisGamman))
        
                # Define Lgas for input ranges
                Lgas[idx, ifreq] = interpLosses(Rkm, Rs, angamman)
        else:
            # Compute cumulative loss for layers
            thisGamman = np.hstack((0, interp1(midpointLayers[idxNotNaN], gamman[idxNotNaN], hlayers)))
            thisGamman[isnan(thisGamman)] = 0
            angamman = cumsum(bsxfun("times", an, thisGamman))
    
            # Define Lgas for input ranges
            Lgas[idx, 0] = interpLosses(Rkm, Rs, angamman)
    
    # Calculate lens loss (optional)
    # if nargout == 2:
    Llens = lenspl(R, platformHeight, el, VaporDensity = 7.5, ScaleHeight = 2000, LatitudeModel = "Standard", Season = "Summer", AtmosphereMeasurements = None)
    # varargout{1} = repmat(Llens, 1, numF)
    varargout = repmat(Llens, 1, numF)
    
    return Lgas, varargout
    
def tropopl_getRange(R, ie, numEl):
    Rcol = singcolvec(R)
    if numEl == 1:
        idx = list(range(0, numel(R)))
    else:
        idx = ie

    Rkm = singrowvec(Rcol[idx, :] * 1e-3)
    
    return Rkm, idx
    
def interpLosses(Rkm, Rs, angamman):
    Rkm = singrowvec(Rkm)
    Rs = singrowvec(Rs)

    idxG = Rkm > matlabmax(Rs, [], 2)[0]
    
    # Calculate loss for specified range values
    L = zeros(numel(Rkm))
    if np.any(~idxG):
        L[~idxG] = interp1(Rs, angamman, Rkm[~idxG], 'linear')

    if np.any(idxG):
        L[idxG] = angamman(-1)

    return L

def calculatelayers(modelStr, tpwvdenh):
    '''
       This function is for internal use only. It may be removed in the
       future.
    
     calculatelayers    Calculate ITU recommended layers for refraction
       EDGES = calculatelayers(MODELSTR,TPWVDENH) calculates the ITU
       recommended thickness of layers for a sufficient estimation of
       refraction and path attenuation. The thickness of the layers in this
       model increases exponentially from 10 cm at ground level to about 1 km
       at about 100 km.
    
       The input MODELSTR is a string that corresponds to the atmosphere model
       used. If the model string is not 'Custom', an ITU reference atmosphere
      is assumed, which has a minimum altitude of 0 km MSL. If a custom model
       is used, the minimum altitude is set according to the first value in
       the last column of TPWVDENH, which is defined as per the command-line
       help in refractiveidx, tropopl, and lenspl.
    
       The output EDGES are the layer edges in km.
    
       Copyright 2020 The MathWorks, Inc.
    
       Reference 
       [1] International Telecommunication Union (ITU). "Attenuation by
           Atmospheric Gases." Recommendation ITU-R P.676-11, P Series,
           Radiowave Propagation, Sept. 2016.
    '''
    # Set minimum and maximum altitude
    if strcmpi(modelStr, 'custom'):
        minAlt = tpwvdenh[0, 3] * 1e-3 # km
    else:
        minAlt = 0 # km
    
    # Define ITU layers
    idx = np.array(list(range(0, 922)))

    if minAlt == 0:
        edges = cumsum(np.concatenate((np.array([0]), 0.0001 * exp(idx / 100)), axis=0)) # Eqn 21
    else:
        edges = cumsum(np.concatenate((np.arange(minAlt, 0, 0.0001), 0.0001 * exp(idx / 100)), axis=0)) # Eqn 21

    return edges

def gasatt(f, T, P, rho):
    '''
    %   This function is for internal use only. It may be removed in the
    %   future.
    
    % gasatt    Calculates the tropospheric gamma attenuation factor
    %   GAMMA = gasatt(F,T,P,RHO) calculates the tropospheric attenuation
    %   factor GAMMA, which is an M x N matrix where M is the number of
    %   frequencies and N is the number of tropospheric layers. 
    %
    %   The input F is the an M-length vector of frequencies in Hz, T is a 1xN
    %   vector of the tropospheric temperatures in Kelvin for layer n, P is a
    %   1xN vector of the tropospheric pressures in hPa for layer n, and rho is
    %   a 1xN vector of the water vapor densities in g/m^3 for layer n.
    
    %   Reference 
    %   [1] International Telecommunication Union (ITU). "Attenuation by
    %       Atmospheric Gases." Recommendation ITU-R P.676-12, P Series,
    %       Radiowave Propagation, Aug. 2019.
    
    %   Copyright 2020-2023 The MathWorks, Inc.
    '''
    
    # Convert frequency to GHz 
    fGHz = f / 1e9
    
    # Calculate theta
    theta = 300 / T
    
    # Calculate water vapor pressure 
    Pwv = rho * T / 216.7 # hPa 
    
    # Calculate dry air and water vapor coefficients 
    gammao = dryaircoeff(fGHz, theta, P, Pwv)
    gammaw = watvapcoeff(fGHz, theta, P, Pwv)
    
    # Combined specific attenuation
    gamma = gammao + gammaw
    
    return gamma
    
def dryaircoeff(fGHz, theta, P, Pwv):
    '''
       This function is for internal use only. It may be removed in the
       future.
    
       gammao = dryaircoeff(fGHz,theta,P,Pwv) calculates specific attenuation
       for dry air according to ITU-R P.676-12.
    
       The input fGHz is an Mx1 vector of frequencies in GHz, theta is a 1xN
       length vector equal to 300 divided by the temperature in Kelvin, P is a
       1xN vector of the tropospheric pressures in hPa for layer n, and Pwv is
       a 1xN vector of the water vapor pressures in hPa for layer n.
    
       The output gammo is a MxN matrix where M corresponds to the number of
       frequencies and N corresponds to the number of layers. 
    
       Reference 
       [1] International Telecommunication Union (ITU). "Attenuation by
           Atmospheric Gases." Recommendation ITU-R P.676-12, P Series,
           Radiowave Propagation, Aug. 2019.
    
       Copyright 2020-2021 The MathWorks, Inc.
    '''
    
    oxytab = np.array([
        # f0 ,a1 ,a2 ,a3 ,a4 ,a5 ,a6
        [50.474214, 0.975, 9.651, 6.690, 0.0, 2.566, 6.850],
        [50.987745, 2.529, 8.653, 7.170, 0.0, 2.246, 6.800],
        [51.503360, 6.193, 7.709, 7.640, 0.0, 1.947, 6.729],
        [52.021429, 14.320, 6.819, 8.110, 0.0, 1.667, 6.640],
        [52.542418, 31.240, 5.983, 8.580, 0.0, 1.388, 6.526],
        [53.066934, 64.290, 5.201, 9.060, 0.0, 1.349, 6.206],
        [53.595775, 124.600, 4.474, 9.550, 0.0, 2.227, 5.085],
        [54.130025, 227.300, 3.800, 9.960, 0.0, 3.170, 3.750],
        [54.671180, 389.700, 3.182, 10.370, 0.0, 3.558, 2.654],
        [55.221384, 627.100, 2.618, 10.890, 0.0, 2.560, 2.952],
        [55.783815, 945.300, 2.109, 11.340, 0.0, -1.172, 6.135],
        [56.264774, 543.400, 0.014, 17.030, 0.0, 3.525, -0.978],
        [56.363399, 1331.800, 1.654, 11.890, 0.0, -2.378, 6.547],
        [56.968211, 1746.600, 1.255, 12.230, 0.0, -3.545, 6.451],
        [57.612486, 2120.100, 0.910, 12.620, 0.0, -5.416, 6.056],
        [58.323877, 2363.700, 0.621, 12.950, 0.0, -1.932, 0.436],
        [58.446588, 1442.100, 0.083, 14.910, 0.0, 6.768, -1.273],
        [59.164204, 2379.900, 0.387, 13.530, 0.0, -6.561, 2.309],
        [59.590983, 2090.700, 0.207, 14.080, 0.0, 6.957, -0.776],
        [60.306056, 2103.400, 0.207, 14.150, 0.0, -6.395, 0.699],
        [60.434778, 2438.000, 0.386, 13.390, 0.0, 6.342, -2.825],
        [61.150562, 2479.500, 0.621, 12.920, 0.0, 1.014, -0.584],
        [61.800158, 2275.900, 0.910, 12.630, 0.0, 5.014, -6.619],
        [62.411220, 1915.400, 1.255, 12.170, 0.0, 3.029, -6.759],
        [62.486253, 1503.000, 0.083, 15.130, 0.0, -4.499, 0.844],
        [62.997984, 1490.200, 1.654, 11.740, 0.0, 1.856, -6.675],
        [63.568526, 1078.000, 2.108, 11.340, 0.0, 0.658, -6.139],
        [64.127775, 728.700, 2.617, 10.880, 0.0, -3.036, -2.895],
        [64.678910, 461.300, 3.181, 10.380, 0.0, -3.968, -2.590],
        [65.224078, 274.000, 3.800, 9.960, 0.0, -3.528, -3.680],
        [65.764779, 153.000, 4.473, 9.550, 0.0, -2.548, -5.002],
        [66.302096, 80.400, 5.200, 9.060, 0.0, -1.660, -6.091],
        [66.836834, 39.800, 5.982, 8.580, 0.0, -1.680, -6.393],
        [67.369601, 18.560, 6.818, 8.110, 0.0, -1.956, -6.475],
        [67.900868, 8.172, 7.708, 7.640, 0.0, -2.216, -6.545],
        [68.431006, 3.397, 8.652, 7.170, 0.0, -2.492, -6.600],
        [68.960312, 1.334, 9.650, 6.690, 0.0, -2.773, -6.650],
        [118.750334, 940.300, 0.010, 16.640, 0.0, -0.439, 0.079],
        [368.498246, 67.400, 0.048, 16.400, 0.0, 0.000, 0.000],
        [424.763020, 637.700, 0.044, 16.400, 0.0, 0.000, 0.000],
        [487.249273, 237.400, 0.049, 16.000, 0.0, 0.000, 0.000],
        [715.392902, 98.100, 0.145, 16.000, 0.0, 0.000, 0.000],
        [773.839490, 572.300, 0.141, 16.200, 0.0, 0.000, 0.000],
        [834.145546, 183.100, 0.145, 14.700, 0.0, 0.000, 0.000],
    ])
    
    M = numel(fGHz)
    N = numel(theta)
    
    # Eqn 3
    # Si = Line strength
    # Si = a1*10^-7*P*theta^3*exp(a2*(1 - theta))
    Si = bsxfun("times", oxytab[:, 1] * 1e-7, P * theta ** 3) * exp(bsxfun("times", oxytab[:, 2], (1 - theta))) # 3
    
    # Eqn 6a and 6b
    # deltaf = a3*10^-4*(P*theta^(0.8 - a4) + 1.1*Pwv*theta)
    # deltaf = sqrt(deltaf^2 + 2.25*10^-6)
    deltafvecA = bsxfun("times", P, bsxfun("power", theta, 0.8 - oxytab[:, 4]))
    deltafvecB = bsxfun("plus", deltafvecA, 1.1 * Pwv * theta)
    deltafvec = bsxfun("times", oxytab[:, 3] * 1e-4, deltafvecB) # 6a
    deltafvec = sqrt(deltafvec ** 2 + 2.25e-6) # 6b
    
    # Eqn 7
    # delta = (a5 + a6*theta)*10^-4*(P + Pwv)*theta^0.8
    deltavec = bsxfun("times", bsxfun("plus", oxytab[:, 5], bsxfun("times", oxytab[:, 6], theta)) * 1e-4, (P + Pwv) * theta ** 0.8) # 7
    
    # # Eqn 5 
    # Fi = Line factor
    # # Fi = f/fi*((deltaf - delta*(fi - f))/((fi - f)^2 + deltaf^2) + (deltaf -
    # # delta*(fi + f))/((fi + f)^2 + deltaf^2))
    fdiff = bsxfun("minus", oxytab[:, 0], fGHz)#.')
    fsum = bsxfun("plus", oxytab[:, 0], fGHz)#.')
    if np.size(fGHz) > 1:
        sumSiFiO = zeros((M, N))
        for im in range(0, M):
            FiPart = (deltafvec - bsxfun("times", deltavec, fdiff[:, im])) / (bsxfun("plus", fdiff[:, im] ** 2, deltafvec ** 2)) + (deltafvec - bsxfun("times", deltavec, fsum[:, im])) / (bsxfun("plus", fsum[:, im] ** 2, deltafvec ** 2))
            Fi = bsxfun("times", fGHz[im] / oxytab[:, 0], FiPart) 
            sumSiFiO[im, :] = matlabsum(Si * Fi, 1)
    else:
        sumSiFiO = zeros(N)
        FiPart = (deltafvec - bsxfun("times", deltavec, fdiff)) / (bsxfun("plus", fdiff ** 2, deltafvec ** 2)) + (deltafvec - bsxfun("times", deltavec, fsum)) / (bsxfun("plus", fsum ** 2, deltafvec ** 2))
        Fi = bsxfun("times", fGHz / oxytab[:, 0], FiPart)
        sumSiFiO = matlabsum(Si * Fi, 1)

    # Calculate specific attenuation for oxygen
    d = 5.6e-4 * (P + Pwv) * theta ** 0.8 # Eqn 9
    # NDdp = fGHz.*P.*theta.^2.*(6.14e-5./(d.*(1+(fGHz./d).^2)) + ...
    #     1.4e-12.*P.*theta.^1.5./(1 + 1.9e-5.*fGHz.^1.5)); % Eqn 8
    NDdp = bsxfun("times", fGHz, P * theta ** 2) * (6.14e-5 / bsxfun("times", d, 1 + (bsxfun("rdivide", fGHz, d)) ** 2) + bsxfun("rdivide", 1.4e-12 * P * theta ** 1.5, (1 + 1.9e-5 * fGHz ** 1.5))) # Eqn 8
    gammao = 0.1820 * bsxfun("times", fGHz, (sumSiFiO + NDdp)) # Specific attenuation for oxygen

    return gammao

def watvapcoeff(fGHz, theta, P, Pwv):
    '''
    %   This function is for internal use only. It may be removed in the
    %   future.
    
    %   gammaw = watvapcoeff(fGHz,theta,P,Pwv) calculates specific attenuation
    %   for water according to ITU-R P.676-12.
    %
    %   The input fGHz is an Mx1 vector of frequencies in GHz, theta is a 1xN
    %   length vector equal to 300 divided by the temperature in Kelvin, P is a
    %   1xN vector of the tropospheric pressures in hPa for layer n, and Pwv is
    %   a 1xN vector of the water vapor pressures in hPa for layer n.
    %
    %   The output gammw is a MxN matrix where M corresponds to the number of
    %   frequencies and N corresponds to the number of layers.
    
    %   Reference 
    %   [1] International Telecommunication Union (ITU). "Attenuation by
    %       Atmospheric Gases." Recommendation ITU-R P.676-12, P Series,
    %       Radiowave Propagation, Aug. 2019.
    
    %   Copyright 2020-2021 The MathWorks, Inc.
    '''    
    watvaptab = np.array([
        # f0, b1, b2, b3, b4, b5, b6
        [22.235080, 0.1079, 2.144, 26.38, 0.76, 5.087, 1.00],
        [67.803960, 0.0011, 8.732, 28.58, 0.69, 4.930, 0.82],
        [119.995940, 0.0007, 8.353, 29.48, 0.70, 4.780, 0.79],
        [183.310087, 2.273, 0.668, 29.06, 0.77, 5.022, 0.85],
        [321.225630, 0.0470, 6.179, 24.04, 0.67, 4.398, 0.54],
        [325.152888, 1.514, 1.541, 28.23, 0.64, 4.893, 0.74],
        [336.227764, 0.0010, 9.825, 26.93, 0.69, 4.740, 0.61],
        [380.197353, 11.67, 1.048, 28.11, 0.54, 5.063, 0.89],
        [390.134508, 0.0045, 7.347, 21.52, 0.63, 4.810, 0.55],
        [437.346667, 0.0632, 5.048, 18.45, 0.60, 4.230, 0.48],
        [439.150807, 0.9098, 3.595, 20.07, 0.63, 4.483, 0.52],
        [443.018343, 0.1920, 5.048, 15.55, 0.60, 5.083, 0.50],
        [448.001085, 10.41, 1.405, 25.64, 0.66, 5.028, 0.67],
        [470.888999, 0.3254, 3.597, 21.34, 0.66, 4.506, 0.65],
        [474.689092, 1.260, 2.379, 23.20, 0.65, 4.804, 0.64],
        [488.490108, 0.2529, 2.852, 25.86, 0.69, 5.201, 0.72],
        [503.568532, 0.0372, 6.731, 16.12, 0.61, 3.980, 0.43],
        [504.482692, 0.0124, 6.731, 16.12, 0.61, 4.010, 0.45],
        [547.676440, 0.9785, 0.158, 26.00, 0.70, 4.500, 1.00],
        [552.020960, 0.1840, 0.158, 26.00, 0.70, 4.500, 1.00],
        [556.935985, 497.0, 0.159, 30.86, 0.69, 4.552, 1.00],
        [620.700807, 5.015, 2.391, 24.38, 0.71, 4.856, 0.68],
        [645.766085, 0.0067, 8.633, 18.00, 0.60, 4.000, 0.50],
        [658.005280, 0.2732, 7.816, 32.10, 0.69, 4.140, 1.00],
        [752.033113, 243.4, 0.396, 30.86, 0.68, 4.352, 0.84],
        [841.051732, 0.0134, 8.177, 15.90, 0.33, 5.760, 0.45],
        [859.965698, 0.1325, 8.055, 30.60, 0.68, 4.090, 0.84],
        [899.303175, 0.0547, 7.914, 29.85, 0.68, 4.530, 0.90],
        [902.611085, 0.0386, 8.429, 28.65, 0.70, 5.100, 0.95],
        [906.205957, 0.1836, 5.110, 24.08, 0.70, 4.700, 0.53],
        [916.171582, 8.400, 1.441, 26.73, 0.70, 5.150, 0.78],
        [923.112692, 0.0079, 10.293, 29.00, 0.70, 5.000, 0.80],
        [970.315022, 9.009, 1.919, 25.50, 0.64, 4.940, 0.67],
        [987.926764, 134.6, 0.257, 29.85, 0.68, 4.550, 0.90],
        [1780.000000, 17506, 0.952, 196.3, 2.00, 24.15, 5.00]
    ])
    
    M = numel(fGHz)
    N = numel(theta)
    
    # Eqn 3
    # Si = b1*10^-1*Pwv*theta^3.5*exp(b2*(1 - theta))
    Si = bsxfun("times", watvaptab[:, 1] * 1e-1, Pwv * theta ** 3.5) * exp(bsxfun("times", watvaptab[:, 2], (1 - theta))) # eq(3)
    
    # Eqn 6a and 6b
    # deltaf = b3*10^-4*(P*theta^b4 + b5*Pwv*theta^b6)
    # deltaf = 0.535*deltaf + sqrt(0.217*deltaf^2 + 2.1316*10^-12*fi^2/theta)
    deltafvecA = bsxfun("times", P, bsxfun("power", theta, watvaptab[:, 4]))
    deltafvecB = bsxfun("times", watvaptab[:, 5], bsxfun("times", Pwv, bsxfun("power", theta, watvaptab[:, 6])))
    deltafvec = bsxfun("times", watvaptab[:, 3] * 1e-4, bsxfun("plus", deltafvecA, deltafvecB)) # 6a
    deltafvec = 0.535 * deltafvec + sqrt(0.217 * deltafvec ** 2 + 2.1316e-12 * bsxfun("rdivide", watvaptab[:, 0] ** 2, theta)) # 6b
    
    # Eqn 7
    # delta = 0
    deltavec = zeros(size(deltafvec)) # eq(7)
    
    # Eqn 5 
    # Fi = f/fi*((deltaf - delta*(fi - f))/((fi - f)^2 + deltaf^2) + (deltaf -
    # delta*(fi + f))/((fi + f)^2 + deltaf^2))
    fdiff = bsxfun("minus", watvaptab[:, 0], fGHz)#.');
    fsum = bsxfun("plus", watvaptab[:, 0], fGHz)#.');
    sumSiFiW = zeros((M, N))
    if np.size(fGHz) > 1:
        for im in range(0, M):
            FiPart = (deltafvec - bsxfun("times", deltavec, fdiff[:, im])) / (bsxfun("plus", fdiff[:, im] ** 2, deltafvec ** 2)) + (deltafvec - bsxfun("times", deltavec, fsum[:, im])) / (bsxfun("plus", fsum[:, im] ** 2, deltafvec ** 2))
            Fi = bsxfun("times", fGHz[im] / watvaptab[:, 0], FiPart) # eq(5)
            sumSiFiW[im, :] = matlabsum(Si * Fi, 1)
    else:
        FiPart = (deltafvec - bsxfun("times", deltavec, fdiff)) / (bsxfun("plus", fdiff ** 2, deltafvec ** 2)) + (deltafvec - bsxfun("times", deltavec, fsum)) / (bsxfun("plus", fsum ** 2, deltafvec ** 2))
        Fi = bsxfun("times", fGHz / watvaptab[:, 0], FiPart) # eq(5)
        sumSiFiW = matlabsum(Si * Fi, 1)
    
    # Calculate specific attenuation for water
    gammaw = 0.1820 * bsxfun("times", fGHz, sumSiFiW) # Specific attenuation for water

    return gammaw

def lenspl(R, platformHeight, el, VaporDensity = 7.5, ScaleHeight = 2000, LatitudeModel = "Standard", Season = "Summer", AtmosphereMeasurements = None):
    '''
    % lenspl    Calculates the loss due to the tropospheric lens effect 
    %   L = lenspl(R,H,EL) calculates the one-way loss due to the tropospheric
    %   lens effect using the International Telecommunication Union (ITU)
    %   standard atmospheric model known as the Mean Annual Global Reference
    %   Atmosphere (MAGRA), which approximates the U.S. Standard Atmosphere
    %   1976 with insignificant relative error. The variation in refractivity
    %   versus altitude makes the atmosphere act like a lens with loss
    %   independent of frequency. Rays leaving an antenna are refracted in the
    %   troposphere and the energy radiated within some angular extent is
    %   distributed over a slightly greater angular sector, thereby reducing
    %   the energy density relative to propagation in a vacuum.
    %
    %   H is a scalar indicating the mean sea level (MSL) altitude of the radar
    %   platform in m, EL specifies the scalar or M-length vector of initial
    %   elevation angles (in deg) at the radar platform, and R is an M-length
    %   vector indicating the slant range in m. The output L is an M-length
    %   column vector of the loss in dB.
    %
    %   L = lenspl(...,'WaterVaporDensity',RHO0) specifies the standard
    %   ground-level water vapor density RHO0 as a nonnegative scalar in g/m^3.
    %   Applicable only for the default standard model (Mean Annual Global
    %   Reference Atmosphere). Defaults to 7.5 g/m^3.
    %
    %   L = lenspl(...,'ScaleHeight',H0) specifies the scale height H0
    %   (altitude above MSL) as a nonnegative scalar in m. Applicable only for the
    %   default standard model (Mean Annual Global Reference Atmosphere).
    %   Defaults to 2e3 m. For a dry atmosphere, set H0 to 6e3 m.
    %
    %   L = lenspl(...,'LatitudeModel',MODEL) specifies the reference
    %   latitude model used. Defaults to 'Standard'. Applicable string inputs
    %   are as follows:
    %      - 'Standard' 
    %           This model is the Mean Annual Global Reference Atmosphere
    %           (MAGRA) that reflects the mean annual temperature and pressure
    %           averaged across the world. Default model.
    %      - 'Low' 
    %           This model is for low latitudes less than 22 deg, where there
    %           exists little seasonal variation.
    %      - 'Mid' 
    %           This model is for mid-latitudes between 22 deg and 45 deg with
    %           seasonal profiles for 'Summer' and 'Winter', which can be
    %           specified using the 'Season' name-value pair as discussed more
    %           below. Defaults to 'Summer' profile.
    %      - 'High' 
    %           This model is for high-latitudes greater than 45 deg with
    %           seasonal profiles for 'Summer' and 'Winter', which can be
    %           specified using the 'Season' name-value pair as discussed more
    %           below. Defaults to 'Summer' profile.
    %
    %   L = lenspl(...,'Season',SEASON) specifies the season as a string. Valid
    %   only for the 'Mid' and 'High' latitude models; other models will ignore
    %   this input. Applicable inputs are 'Summer' or 'Winter'. Defaults to
    %   'Summer'.
    %
    %   L = lenspl(...,'AtmosphereMeasurements',TPWVDENH) specifies custom
    %   atmospheric measurements for the calculation of the refractive index.
    %   TPWVDENH is an N x 4 matrix, where N corresponds to the number of
    %   altitude measurements. The first column is the atmospheric temperature
    %   in K, the second column is the atmospheric dry air pressure in hPa, the
    %   third column is the water vapor density in g/m^3, and the fourth column
    %   is the MSL altitude of the measurements in m. When a custom model is
    %   provided, all other name-value pair options are ignored and the output
    %   refractive index is applicable for the input height H.
    %
    %   Lens loss is not due to a dissipation of energy and does not contribute
    %   to the system noise temperature. Lens loss is a function of range and
    %   elevation angle, and the loss varies slowly at long ranges. Lens loss
    %   can be significant at long ranges and shallow propagation angles. This
    %   function only calculates lens loss for positive elevation angles, as
    %   per the original paper. Negative elevation angles return loss equal to
    %   0. 
    %
    %   Examples:
    %
    %   % Example 1:
    %   %   Calculate the one-way lens loss for a platform height of 10 m with 
    %   %   a grazing angle of 0.5 deg at a range of 200 km. 
    %   L = lenspl(200e3,10,0.5) 
    %   
    %   % Example 2: 
    %   %   Calculate the one-way lens loss for a platform height of 10 m with 
    %   %   a grazing angle of 0.5 deg at a range of 200 km using the 
    %   %   low-latitude refractivity model.
    %   L = lenspl(200e3,10,0.5,'LatitudeModel','Low') 
    %   
    %   % Example 3:
    %   %   Calculate and plot the two-way lens loss curve at an elevation 
    %   %   angle of 0.03 deg. 
    %   R  = (100:5000).*1e3;  % m
    %   h  = 0;                % m
    %   el = 0.03;             % deg
    %   L  = 2*lenspl(R,h,el); % Factor of 2 for two-way propagation
    %   plot(R.*1e-3,L)
    %   grid on; 
    %   xlabel('Range (km)')
    %   ylabel('Loss (dB)')
    %   title('Two-Way Lens Loss')
    %
    %   See also atmositu, refractiveidx, tropopl, gaspl, effearthradius.
    
    %   Copyright 2020-2022 The MathWorks, Inc.
    
    %   Reference
    %   [1] Weil, T. A. "Atmospheric Lens Effect: Another Loss for the Radar
    %       Range Equation." IEEE Transactions on Aerospace and Electronic
    %       Systems, Vol. AES-9, No. 1, Jan. 1973.
    %   [2] Meikle, H. Modern Radar Systems. 2nd edition. Norwood, MA: Artech
    %       House, 2008.
    %   [3] Barton, David K. Radar Equations for Modern Radar. 1st edition.
    %       Norwood, MA: Artech House, 2013.
    %   [4] International Telecommunication Union (ITU). "Attenuation by
    %       Atmospheric Gases." Recommendation ITU-R P.676-11, P Series,
    %       Radiowave Propagation, Sept. 2016.
    '''
    
    hgeomKm = platformHeight * 1e-3 # Geometric height (km)
    T, P, wvden = stdatm(hgeomKm, VaporDensity, ScaleHeight)

    if AtmosphereMeasurements is None:
        tpwvdenh = [singcolvec(T), singcolvec(P), singcolvec(wvden), singcolvec(platformHeight)]
    else:
        tpwvdenh = AtmosphereMeasurements
    modelStr = getModelStr(LatitudeModel, Season)
    
    # Lens effect
    #   There are a multitude of ways in which the loss can be calculated.
    #   One way is to express the loss as the derivative of the ray's final
    #   elevation angle with respect to its value as it leaves the antenna,
    #   which is the approach taken here.
    
    # Earth radius
    Re = Rearth * 1e-3 # km
    
    # Calculate layers
    platformHeightKm = platformHeight * 1e-3
    edgesITU0 = calculatelayers(modelStr, tpwvdenh)
    edgesITU = np.hstack((edgesITU0, np.arange(edgesITU0[-1] + 1, int(1e4)))) # Add additional layers for long-range calculations
    midpointLayers = 500 * edgesITU[0:-1] + 500 * edgesITU[1:] # m
    
    # Check for ducting conditions
    ridx, N = refractiveidx(midpointLayers, VaporDensity = VaporDensity, ScaleHeight = ScaleHeight, LatitudeModel = LatitudeModel, Season = Season, AtmosphereMeasurements = AtmosphereMeasurements)
    
    refGrad = diff(N) / diff(midpointLayers)
    cond = any(refGrad <= -157) # From ITU-R P.676-11
    if cond:
        raise ValueError("Lenspl: refgrad below -157")
    
    # Check for subrefractive conditions
    cond = any(refGrad > 0)
    if cond:
        raise ValueError("Lenspl: refgrad cannot be positive")
    
    # This function is intended for shallow grazing angles only. Below about 10
    # degrees, the losses seen are insignificant. 
    L = zeros((numel(R), 1))
    numEl = numel(el)

    if isscalar(el):
        el = [el]
    
    for ie in range(0, numEl):
        if (el[ie] >= 0) and (el[ie] <= 10): # Elevation must be 0 <= el <= 10
            # Below 0.03, numerical errors are present and solutions diverge
            # greatly from expectations
            padVal = 1e-3 * pi / 180 # rad
            if el[ie] < 0.03:
                el[ie] = 0.03

            theta0 = el[ie] * pi / 180 # rad
    
            # Define ideal (no spreading)
            A = 2 * padVal
    
            # Calculate propagation for theta0
            edges0, an0, alphan0, _, _ = propagatePosEl(platformHeightKm, theta0, edgesITU, Re, ridx)

            # Calculate slant range
            Rs, _, _, _ = refractiongeometry(theta0, an0, alphan0, edges0, Re)

            # Diverge +- rays around theta0
            # Rays +/-
            theta0Plus = theta0 + padVal # rad
            theta0Minus = theta0 - padVal # rad
    
            # Calculate propagation for +- rays
            edgesPlus, anPlus, alphanPlus, _,  _ = propagatePosEl(platformHeightKm, theta0Plus, edgesITU, Re, ridx)
            edgesMinus, anMinus, alphanMinus, _, _ = propagatePosEl(platformHeightKm, theta0Minus, edgesITU, Re, ridx)
            
            # Calculate refraction geometry and broadening for +/- rays
            _, thetatPlusGang, _, _,= refractiongeometry(theta0Plus, anPlus, alphanPlus, edgesPlus, Re)
            _, thetatMinusGang, _, _, = refractiongeometry(theta0Minus, anMinus, alphanMinus, edgesMinus, Re)
    
            # Calculate growth in +- rays
            B = abs(thetatPlusGang - thetatMinusGang)
    
            # Calculate
            Rkm, idx = lenspl_getRange(R, ie, numEl)
    
            # Calculate loss curve
            Lcurve = 10 * log10(B / A)
            L[idx] = lenspl_interpLosses(Rkm, Rs, Lcurve).reshape(numel(R), 1)

    return L
    
def lenspl_getRange(R, ie, numEl):
    Rcol = singcolvec(R)
    if numEl == 1:
        idx = list(range(0, numel(R)))
    else:
        idx = ie

    Rkm = Rcol[idx, :] * 1e-3#).';
    return Rkm, idx
    
def lenspl_interpLosses(Rkm, Rs, Lcurve):
    # Calculate loss for specified range values
    Rkm = singrowvec(Rkm)
    L = Nan(numel(Rkm))
    idxLessThan = Rkm <= Rs[-1]

    if any(idxLessThan):
        for idx in idxLessThan:
            L[idx] = interp1(Rs, Lcurve, Rkm[idx], method = 'linear', extrap = True)#.'

    if any(~idxLessThan): # Any greater than max range? Assume maximum loss.
        for idx in ~idxLessThan:
            L[idx] = Lcurve[-1]

    idxLessThan0 = L < 0
    L[idxLessThan0] = 0

    return L

def propagatePosEl(hKm, elRad, edgesITU, Re, ridx):
    '''
    %   This function is for internal use only. It may be removed in the
    %   future.
    
    % propagatePosEl    Propagate positive elevation angle
    %   [EDGES,AN,ALPHAN,BETAN,HLAYERS] =
    %   propagatePosEl(HKM,ELRAD,EDGESITU,RE,RIDX) calculates the propagation
    %   path assuming a positive elevation angle.
    %
    %   HKM is the platform height in km, ELRAD is the positive elevation angle
    %   in radians, EDGESITU are the ITU atmospheric layer edges (or custom
    %   measurements) in km, RE is the Earth radius in km, and RIDX are the
    %   refractive indices corresponding to the midpoints of EDGESITU. 
    %
    %   EDGES are the updated layers based on platform height in km, AN is the
    %   path length through the nth layer in km, ALPHAN is the entry angle in
    %   rad, and BETAN is the exit angle in rad (for testing purposes), and
    %   HLAYERS are the atmospheric layer midpoints used for the refractivity
    %   calculation in m.
    
    %   Copyright 2020-2022 The MathWorks, Inc.
    '''    
    edges = updatelayers(hKm, edgesITU) # km
    hlayers = 500 * edges[0:-1] + 500 * edges[1:] # m
    edgesMidpoint = 500 * edgesITU[0:-1] + 500 * edgesITU[1:] # m
    ridxn = interp1(edgesMidpoint, ridx, hlayers)
    ridxn[isnan(ridxn)] = 1
    
    # Refract ray
    an, alphan, betan, _ = refractray(elRad, ridxn, edges, Re)
    
    return edges, an, alphan, betan, hlayers

def updatelayers(platformHeight, edgesIn):
    '''
    %   This function is for internal use only. It may be removed in the
    %   future.
    
    % updatelayers    Updates ITU layers based on platform height
    %   [EDGES] = updatelayers(H,EDGESin) updates the refraction layers
    %   based on the platform altitude.
    %
    %   H is the MSL altitude of the platform in km and EDGES are the layer
    %   edges in km.
    %   
    %   EDGES is a subset of the corresponding input in km. 
    
    %   Copyright 2020 The MathWorks, Inc.
    '''
    edgesIn = numpyify(edgesIn)
    if platformHeight == 0:
        # Calculation is for all layers
        edges = edgesIn
    else: # platformHeight < maxAlt
        # Calculation is a subset of layers
        nlayers = numel(edgesIn) - 1
        idxAlt = interp1(edgesIn[0:len(edgesIn) - 1], list(range(0, nlayers)), platformHeight, 'linear') # Potentially non-integer
        idxAltC = int(ceil(idxAlt)) # Integer index ceil-ed
        if mod(idxAlt, 1) == 0:
            edges = edgesIn[idxAlt:]
        else:
            edges = [platformHeight, *edgesIn[idxAltC:]]

    return edges

def refractray(el, ridxn, edgesIn, Re):
    tol = 1e-3  # Tolerances for isnearlyequal check
    zenithcase = False
    if abs(el) > 0:
        tf, el = isnearlyequal(el, copysign(pi / 2, el), tol)
        zenithcase = tf

    # Initialize layer-by-layer calculation
    deltan = [abs(edgesIn[i+1] - edgesIn[i]) for i in range(len(edgesIn) - 1)]
    nlayers = len(deltan)

    if zenithcase:
        alphan = [0.0] * nlayers
        betan = [0.0] * nlayers
        an = deltan  # km, Equation 17
        edges = edgesIn
    elif el >= 0:
        # Initialize
        alphan = [0.0] * nlayers
        betan = [0.0] * nlayers
        rn = [edgesIn[i] + Re for i in range(nlayers)]  # km
        betan[0] = pi / 2 - abs(el)  # rad

        # Calculate angles
        for in_ in range(nlayers):
            alphan[in_] = asin(rn[in_] / (rn[in_] + deltan[in_]) * sin(betan[in_]))  # Eqn 18b
            if in_ < nlayers - 1:
                betan[in_ + 1] = asin(ridxn[in_] / ridxn[in_ + 1] * sin(alphan[in_]))  # Eqn 19a

        # Calculate path length within layers
        an = [
            -rn[i] * cos(betan[i]) + sqrt((rn[i] ** 2 * cos(betan[i]) ** 2 + 2 * rn[i] * deltan[i]) + deltan[i] ** 2)
            for i in range(nlayers)
        ]  # km, Equation 17
        edges = edgesIn
    else:  # el < 0
        # Initialize
        alphanFP = [0.0] * nlayers
        betanFP = [0.0] * nlayers
        idxEndFP = nlayers  # Index of the end of the first path

        # Calculate angles; derived from geometrical theory
        rn = [e + Re for e in edgesIn]  # km
        alphanFP[0] = pi / 2 - abs(el)  # rad
        for in_ in range(nlayers):
            val = rn[in_] / rn[in_ + 1] * sin(alphanFP[in_])
            if val > 1:
                idxEndFP = in_  # Grazing height has been found
                betanFP[in_] = pi / 2
                break
            betanFP[in_] = asin(val)  # rad
            if in_ != nlayers - 1:
                alphanFP[in_ + 1] = asin(ridxn[in_] / ridxn[in_ + 1] * sin(betanFP[in_]))  # rad

        # Calculate path length within layers
        anFP = [
            -rn[i + 1] * cos(betanFP[i]) + sqrt((rn[i + 1] ** 2 * cos(betanFP[i]) ** 2 + 2 * rn[i + 1] * deltan[i]) + deltan[i] ** 2)
            for i in range(nlayers)
        ]  # km, Equation 17

        # Set values for path
        an = anFP[:idxEndFP + 1]
        alphan = alphanFP[:idxEndFP + 1]
        betan = betanFP[:idxEndFP + 1]
        edges = edgesIn[:idxEndFP + 2]

    return an, alphan, betan, edges

def isnearlyequal(a, b, tol):
    a = np.array(a)
    b = np.array(b)

    absA = np.abs(b)
    absB = np.abs(a)
    diffAB = np.abs(a - b)

    tf = np.ones_like(a, dtype=bool)

    idx2 = (a == 0) | (b == 0) | ((absA + absB) < np.finfo(float).tiny)
    tf[idx2] = diffAB[idx2] < tol

    idx3 = ~(a == b) & ~idx2
    tf[idx3] = diffAB[idx3] / np.minimum((absA[idx3] + absB[idx3]), np.finfo(float).max) < tol

    val = np.where(tf, b, a)

    return tf, val

# # Helper function to check if values are nearly equal
# def isnearlyequal(a, b, tol):
#     return abs(a - b) <= tol, b

def propagateNegEl(hKm, elRad, edgesITU, Re, ridx):
    '''
    %   This function is for internal use only. It may be removed in the
    %   future.
    
    % propagateNegEl    Propagate negative elevation angle
    %   [EDGES,AN,ALPHAN,BETAN,HLAYERS] =
    %   propagatePosEl(HKM,ELRAD,EDGESITU,RE,RIDX) calculates the propagation
    %   path assuming a negative elevation angle.
    %
    %   HKM is the platform height in km, ELRAD is the negative elevation angle
    %   in radians, EDGESITU are the ITU atmospheric layer edges (or custom
    %   measurements) in km, RE is the Earth radius in km, and RIDX are the
    %   refractive indices corresponding to the midpoints of EDGESITU. 
    %
    %   EDGES are the updated layers based on platform height in km, AN is the
    %   path length through the nth layer in km, ALPHAN is the entry angle in
    %   rad, and BETAN is the exit angle in rad (for testing purposes), and
    %   HLAYERS are the atmospheric layer midpoints used for the refractivity
    %   calculation in m.
    
    %   Copyright 2020-2022 The MathWorks, Inc.
        
    % In some cases of negative elevation (small, negative elevation
    % angles), there exists two propagation paths. In such a case the wave
    % will propagate along the first path until a minimum grazing height is
    % reached. Then, the wave will propagate in a second path between the
    % minimum grazing height and space. Calculate both paths, and combine
    % to form total path.
    '''
    
    ### Path 1
    # Descending path
    if np.ndim(edgesITU) == 1:
        edgesITU = edgesITU[np.newaxis, :]
        
    edgesD = updatelayers(hKm, fliplr(edgesITU))
    hlayersD = 500 * edgesD[0:-1] + 500 * edgesD[1:] # m
    edgesMidpoint = 500 * np.ravel(edgesITU)[0:-1] + 500 * np.ravel(edgesITU)[1:] # m
    ridxnD = interp1(edgesMidpoint, ridx, hlayersD)
    ridxnD[isnan(ridxnD)] = 1
    
    # Refract ray
    anFP, alphanFP, betanFP, edgesFP = refractray(elRad,ridxnD,edgesD,Re);
    
    ### Path 2
    numLayersFP = numel(anFP)
    if numLayersFP < numel(hlayersD):
        # There is a second path. The first path reached the minimum
        # grazing height.
        
        # Ascending path
        edgesA = updatelayers(edgesFP[-1], np.ravel(edgesITU))
        hlayersA = 500 * edgesA[0:-1] + 500 * edgesA[1:] # m
        ridxnA = interp1(edgesMidpoint, ridx, hlayersA)
        ridxnA[isnan(ridxnA)] = 1
        
        # Refract ray
        elRadSP = 0
        anSP, alphanSP, betanSP, _ = refractray(elRadSP, ridxnA, edgesA, Re)
        
        # Combine
        edges = np.hstack((edgesD[0:numLayersFP], edgesA))
        an = np.hstack((anFP, anSP))
        alphan = np.hstack((alphanFP, alphanSP)) 
        betan = np.hstack((betanFP, betanSP))
        hlayers = np.hstack((hlayersD[0:numLayersFP], hlayersA))
    else:
        # The elevation angle was such that there is not a second path. The
        # wave hit the ground.
        edges = edgesD
        an = anFP
        alphan = alphanFP
        betan = betanFP
        hlayers = hlayersD
    
    return edges, an, alphan, betan, hlayers

def refractiongeometry(el, an, alphan, edges, Re):
    # Check if close to zenith case 
    tol = sqrt(eps(1))  # Tolerances for isnearlyequal check

    zenithcase = False
    if abs(el) > 0:
        # tf = np.isclose(el, np.sign(el)*np.pi/2, atol=tol)
        # if tf:
        #     el = np.sign(el)*np.pi/2
        tf, el = isnearlyequal(el, sign(el) * pi / 2, tol)
        zenithcase = tf

    # Calculate radius from center of earth to beginning of layer n
    if numel(edges) > 1:
        if el < 0:
            rn = edges[1:] + Re  # km
        else:
            rn = edges[:-1] + Re  # km
    else:
        if el < 0:
            rn = edges + Re  # km
        else:
            rn = edges + Re  # km

    if zenithcase:
        # Zenith case
        Rs = cumsum(an)
        thetaE = np.zeros_like(Rs)
        gAng = pi / 2 * np.ones_like(Rs)
        dAng = gAng
    else:
        # Calculate Earth angle
        thetaEn = an * sin(alphan) / rn  # Incremental earth angle changes from rn to rn+1 (rad)
        thetaE = cumsum(thetaEn)  # Earth angle from r1 to rn (rad)

        # Calculation radii 
        r1 = edges[0] + Re  # km
        rNext = edges[1:] + Re  # km

        # Calculate slant range as seen at the platform
        # Note: Limit minimum value to eps to handle numerical issues with small thetaE.
        RsVal = r1 ** 2 + rNext ** 2 - 2 * r1 * rNext * cos(thetaE)
        numRs = numel(RsVal)
        RsVal, _ = matlabmax(np.hstack((singcolvec(RsVal), singcolvec(eps(1) * ones(numRs)))), [], 2)
        Rs = sqrt(singrowvec(RsVal))  # km

        # Configure heights for grazing/depression angle calculation
        if el < 0:
            hsNeg = edges[1:]
            haNeg = edges[0]
            dh = abs(haNeg - hsNeg)
            valG = dh / Rs * (1 + dh / (2 * (Re + hsNeg))) - Rs / (2 * (Re + hsNeg))
            valD = dh / Rs * (1 - dh / (2 * (Re + haNeg))) + Rs / (2 * (Re + haNeg))
        else:
            hsPos = edges[0]
            haPos = edges[1:]
            dh = abs(haPos - hsPos)
            valG = dh / Rs * (1 + dh / (2 * (Re + hsPos))) - Rs / (2 * (Re + hsPos))
            valD = dh / Rs * (1 - dh / (2 * (Re + haPos))) + Rs / (2 * (Re + haPos))

        # Calculate grazing angle changes as seen at the platform
        idx = abs(valG) > 1
        valG[idx] = sign(valG[idx])  # Happens rarely; due to numerical instabilities early on in the iterative calculation
        gAng = asin(valG)

        # Calculate depression angle changes as seen at the platform
        idx = abs(valD) > 1
        valD[idx] = sign(valD[idx])  # Happens rarely; due to numerical instabilities early on in the iterative calculation
        dAng = asin(valD)
        
    return Rs, gAng, dAng, thetaE

def probgrid(p1, p2, n = 100):
    '''
    %probgrid Non-uniformly spaced probabilities
    %   P = probgrid(P1,P2) returns a non-uniformly spaced grid of 100
    %   probabilities between P1 and P2 that correspond to the values of the
    %   normal cumulative distribution function (CDF) evaluated over a set of
    %   points uniformly spaced in the domain of the normal distribution. P1
    %   and P2 must be scalars from a unit interval [0,1] such that P1<P2. The
    %   output P is a row vector.
    %
    %   P = probgrid(P1,P2,N) generates N points between P1 and P2.
    %
    %   % Example:
    %   %   Evaluate the standard normal cumulative distribution function (CDF)
    %   %   on a 10-point grid between 0.2 and 0.95.
    %
    %   pmin = 0.2;
    %   pmax = 0.95;
    %   N = 10;
    %   pd = probgrid(pmin,pmax,N);
    %
    %   % Overlay vector of probabilities on standard normal CDF
    %   x = -3:0.01:3;
    %   sncdf = (1+erf(x/sqrt(2)))/2; % CDF of the standard normal distribution
    %   xd = sqrt(2)*erfinv(2*pd-1);  % Probit function
    %
    %   % Plot
    %   hold on
    %   plot(x,sncdf)
    %   plot(xd,pd,'o')
    %   ylabel('Probability')
    %   legend({'Standard Normal CDF','Probability vector'}, ...
    %       'Location','Northwest')
    %   xticks(xd)
    %   xtickangle(40)
    %   yticks(round(100*pd)/100)
    %   grid on
    %
    %   See also rocinterp, detectability
    
    %   Copyright 2020 The MathWorks, Inc.
    '''
    
    # Step 1: Inverse CDF of the standard normal distribution (a.k.a. "probit"
    # function) used to calculate a uniform grid between xmin and xmax
    xmin = sqrt(2) * erfinv(2 * p1 - 1)
    xmax = sqrt(2) * erfinv(2 * p2 - 1)
    xd = linspace(xmin, xmax, n)
    
    # Step 2: Calculate non-uniform p vector by evaluating the (Standard)
    # Normal CDF on a uniform grid between pmin and pmax
    p = (1 + erf(xd / sqrt(2))) / 2
    p[0] = p1
    p[-1] = p2
    
    return p

def rocinterp(X, Y, Xq, type_str):
    '''
    % rocinterp ROC curve interpolation
    %   PDQ = rocinterp(SNR,PD,SNRQ,'snr-pd') returns the probability of
    %   detection, PDQ, computed by interpolating a probability of detection
    %   vs. signal-to-noise ratio receiver operating characteristic (ROC)
    %   curve. The linear interpolation is performed after transforming the
    %   probability of detection axis of the ROC curve using the normal
    %   probability scale.
    %
    %   The input SNR is a length-J vector of signal-to-noise ratios (in dB)
    %   specifying the horizontal axis of the ROC curve and PD is a length-J
    %   vector or a J-by-K matrix of the corresponding probabilities of
    %   detection. Values in SNR must be unique and PD must be between 0 and 1.
    %   The input SNRQ is a length-N vector specifying the signal-to-noise
    %   ratios at which to interpolate the ROC curve. If PD is a length-J
    %   vector, the output PDQ is a vector of length N. If PD is a J-by-K
    %   matrix, the interpolation is performed along the columns of PD and the
    %   output PDQ is an N-by-K matrix.
    %
    %   SNRQ = rocinterp(PD,SNR,PDQ,'pd-snr') returns the signal-to-noise
    %   ratio, SNRQ, computed by interpolating a probability of detection vs.
    %   signal-to-noise ratio receiver operating characteristic (ROC) curve.
    %   The linear interpolation is performed after transforming the
    %   probability of detection axis of the ROC curve using the normal
    %   probability scale.
    %
    %   The input PD is a length-J vector of probabilities of detection
    %   specifying the vertical axis of the ROC curve and SNR is a length-J
    %   vector or a J-by-K matrix of the corresponding signal-to-noise ratios
    %   (in dB). Values in PD must be unique and must be between 0 and 1. The
    %   input PDQ is a length-N vector specifying the probabilities of
    %   detection at which to interpolate the ROC curve. The values in PDQ must
    %   be between 0 and 1. If SNR is a length-J vector, the output SNRQ is a
    %   vector of length N. If SNR is a J-by-K matrix, the interpolation is
    %   performed along the columns of SNR and the output SNRQ is an N-by-K
    %   matrix.
    %
    %   PDQ = rocinterp(PFA,PD,PFAQ,'pfa-pd') returns the probability of
    %   detection, PDQ, computed by interpolating a probability of detection
    %   vs. probability of false alarm ROC curve. The linear interpolation is
    %   performed after transforming both axes of the ROC curve using a
    %   logarithmic scale.
    %
    %   The input PFA is a length-J vector of probabilities of false alarm
    %   specifying the horizontal axis of the ROC curve and PD is a length-J
    %   vector or a J-by-K matrix of the corresponding probabilities of
    %   detection. Values in PFA must be unique, and both PFA and PD must be
    %   between 0 and 1. The input PFAQ is a length-N vector indicating the
    %   probabilities of false alarm at which to interpolate the ROC curve. The
    %   values in PFAQ must be between 0 and 1. If PD is a length-J vector, the
    %   output PDQ is a vector of length N. If PD is a J-by-K matrix, the
    %   interpolation is performed along the columns of PD and the output PDQ
    %   is an N-by-K matrix.
    %
    %   PFAQ = rocinterp(PD,PFA,PDQ,'pd-pfa') returns the probability of false
    %   alarm, PFAQ, computed by interpolating a probability of detection vs.
    %   probability of false alarm ROC curve. The linear interpolation is
    %   performed after transforming both axes of the ROC curve using a
    %   logarithmic scale.
    %
    %   The input PD is a length-J vector of probabilities of detection
    %   specifying the vertical axis of the ROC curve and PFA is a length-J
    %   vector or a J-by-K matrix of the corresponding probabilities of false
    %   alarm. Values in PD must be unique, and both PD and PFA must be between
    %   0 and 1. The input PDQ is a length-N vector indicating the
    %   probabilities of detection at which to interpolate the ROC curve. The
    %   values in PDQ must be between 0 and 1. If PFA is a length-J vector, the
    %   output PFAQ is a vector of length N. If PFA is a J-by-K matrix, the
    %   interpolation is performed along the columns of PFA and the output PFAQ
    %   is an N-by-K matrix.
    %
    %   % Example:
    %   %   Compute the probability of detection for a Swerling 1 case target
    %   %   given a set of SNR and Pfa values.
    %
    %   pfa = [1e-9 1e-6 1e-3];       % Probability of false alarm
    %   SNR = [13.5 14.5];            % Signal-to-noise ratio (dB)
    %
    %   % Compute the ROC curve
    %   [pd, snr] = rocpfa(pfa,'SignalType','Swerling1');
    %
    %   % Interpolate the ROC curve at the SNR of interest
    %   Pd = rocinterp(snr,pd,SNR,'snr-pd')
    %
    %   See also rocpfa, rocsnr, detectability.
    
    %   Copyright 2020 The MathWorks, Inc.
    
    % When interpolating a ROC curve with varying SNR, map the probability of
    % detection using erfinv, and then perform the linear interpolation. When
    % interpolating a ROC curve with varying Pfa, map probabilities to the log
    % space
    '''
    
    match type_str.lower():
        case 'snr-pd':
            # X is SNR and Y is PD
            Yt = erfinv(2 * Y - 1)
            Yqt = interp1(X, Yt, Xq)
    
            # Map back to probability
            Yq = (1 + erf(Yqt)) / 2
    
        case 'pd-snr':
            # X is PD and Y is SNR
            # Since erfinv maps 1 to Inf and -1 to -Inf;
            X[X > 1-1e-16] = 1-1e-16
            X[X < 1e-16] = 1e-16
    
            Xt = erfinv(2 * X - 1)
            Xqt = erfinv(2 * Xq - 1)
    
            Yq = interp1(Xt, Y, Xqt)
        case _:
            # X is PFA and Y is PD or
            # X is PD and Y is PFA
            Yqt = interp1(log(X), log(Y), log(Xq))
    
            # Map back to probability
            Yq = exp(Yqt)
    
    return Yq

def marcumq(a_in, b_in, m_in = 1):
    '''
    %MARCUMQ Generalized Marcum Q function
    %   MARCUMQ(A,B,M) is the generalized Marcum Q function, defined as
    %
    %      Q_m(a,b) = 1/a^(m-1) * integral from b to inf of
    %                 [x^m * exp(-(x^2+a^2)/2) * I_(m-1)(ax)] dx,
    %
    %   where I_(m-1)() is the modified Bessel function of the first kind, of
    %   order m-1.
    %
    %   a and b must be real and nonnegative.  m must be a positive integer. 
    %   If any of the inputs is a scalar, it is expanded to the size of the 
    %   other inputs.
    %
    %   MARCUMQ(A,B) is the special case for M=1, originally tabulated by
    %   Marcum and often written without the subscript, Q(a,b). 
    %
    %   The calculation in the function is based on the method proposed by
    %   Shnidman, using an absolute error criterion.
    %
    %   % Example:
    %   %   Compute Marcum Q function value of 10 and 20.
    %   
    %   y = marcumq(10,20)
    %
    %   See also BESSELI.
    
    %   Copyright 1996-2020 The MathWorks, Inc.
    
    %   References:
    %     [1] D. A. Shnidman, "The Calculation of the Probability of Detection
    %         and the Generalized Marcum Q-Function", IEEE Transactions on
    %         Information Theory, vol. IT-35, pp. 389-400, March 1989.
    %     [2] J. I. Marcum, "A Statistical Theory of Target Detection
    %         by Pulsed Radar:  Mathematical Appendix", RAND Corporation, Santa
    %         Monica, CA, Research Memorandum RM-753, 1 July 1948. Reprinted in
    %         IRE Transactions on Information Theory, vol. IT-6, pp. 59-267,
    %         April 1960.
    '''
    
    if isscalar(a_in):
        if isscalar(b_in):
            a = repmat(a_in, np.size(m_in))
            b = repmat(b_in, np.size(m_in))
        else:
            b = b_in
            a = repmat(a_in, np.size(b))
            if isscalar(m_in):
                m_in = repmat(m_in, np.size(a))
            else:
                cond = (~isequal(np.size(a), size(m_in)))
                if cond:
                    raise ValueError("marcumq: size of a must equal size of b")

    elif isscalar(b_in):
        a = a_in
        b = repmat(b_in, np.size(a))
        if isscalar(m_in):
            m_in = repmat(m_in, np.size(a))
        else:
            cond = (~isequal(np.size(a), size(m_in)))
            if cond:
                raise ValueError("marcumq: size of a must equal size of b")
    else:
        a = a_in
        b = b_in
        cond = (not isequal(np.size(a), size(b)))
        if cond:
            raise ValueError("marcumq: size of a must equal size of b")
        
        if isscalar(m_in):
            m_in = repmat(m_in, np.size(a));
        else:
            cond = (~isequal(np.size(a), size(m_in)))
            if cond:
                raise ValueError("marcumq: size of a must equal size of b")
    
    if ~isempty(m_in):
        if np.size(m_in) > 1:
            cond = (np.any(m_in < 1) or np.any(m_in / floor(m_in) != 1))
            if cond:
                raise ValueError("marcumq: improper m_in")
        else:
            cond = (m_in < 1) or (m_in / floor(m_in) != 1)
            if cond:
                raise ValueError("marcumq: improper m_in")
    
    # Switch to Shnidman's notation[1].  Note that because the values of a, b,
    # and m_in will affect the choice of power series form and the limits used in
    # the summations evaluated in the P function, we process the vector or
    # matrix one element at a time
    
    Q = Nan(np.size(a))
    for idxx in range(0, numel(a)):
        # Special cases
        if (a[idxx] != inf) and (b[idxx] == 0):
            Q[idxx] = 1
        if np.size(m_in) > 1:
            if all(m_in):
                if (a[idxx] != Inf) and (b[idxx] == Inf):
                    Q[idxx] = 0
                if (a[idxx] == Inf) and (b[idxx] != Inf):
                    Q[idxx] = 1
        else:
            if m_in:
                if (a[idxx] != Inf) and (b[idxx] == Inf):
                    Q[idxx] = 0
                if (a[idxx] == Inf) and (b[idxx] != Inf):
                    Q[idxx] = 1
        
        # General case if the value is still NaN
        # Convert a and b to X and Y using Eq. 7 of [1]
        if np.size(m_in) > 1:
            if isnan(Q[idxx]):
                N = m_in[idxx]
                X = (a[idxx] ** 2 / 2) / N
                Y = b[idxx] ** 2 / 2
                Q[idxx] = P(N, X, Y)
        else:
            if isnan(Q[idxx]):
                N = m_in
                X = (a[idxx] ** 2 / 2) / N
                Y = b[idxx] ** 2 / 2
                Q[idxx] = P(N, X, Y)

    return Q
    
def P(N, X, Y):
    # Computes the probability of detection as given by Shnidman[1].
    
    NX = N * X # Notational convenience
    
    # Evaluate the Chernoff bound
    if Y == 0:
        logCB = 0 # Avoid division by zero
    else:
        Lambda = 1 - N / (2 * Y) - sqrt((N / (2 * Y)) ** 2 + NX / Y) # Eq. 22
        logCB = -Lambda * Y + NX * Lambda / (1 - Lambda) - N * log(1 - Lambda) # Eq. 23
    
    if exp(logCB) < realmin:
        innerArg = 0
        if Y <= NX + N:  
            Q = 1
        else: # if Y > NX + N
            Q = 0
    else:
        innerArg = 0
        if Y < NX + N:   # Compute 1-P(N,X,Y) using form 4 
            Form = 4   # Form 4 corresponds to Eq. 11
        else:            # Compute P(N,X,Y) using form 2 
            Form = 2   # Form 2 corresponds to Eq. 9

        if Form == 2: # wrong
            # Form 2 corresponds to Eq. 9
            # Adjust start of inner summation to avoid underflow
            kmin = findStartOfSummation(NX, 0)
            
            # Adjust start of outer summation to avoid underflow
            mmin = findStartOfSummation(Y, kmin + N)
            
            # Initialize m and k to adjusted values
            m = mmin
            k = m - N
            
            # Initialize the values of innerSum and outerSum
            #
            # If value of k based on N and m is greater than the minimum value
            # computed above, pre-compute start of inner summation.  If not,
            # then we can ignore the first k-1 terms in the inner summation.
            if kmin < k:
                innerSum = 0
                for idx in range(kmin, k):
                    innerArg = expA(NX, idx)
                    if innerArg > realmin:
                        innerSum = innerSum + innerArg
            else:
                innerArg = expA(NX, k)
                innerSum = innerArg
                
            # With m = mmin, ignore the first m-1 terms in the outer summation
            outerArg = expA(Y, m)
            outerSum = outerArg * (1 - innerSum)
            
            # Sum until the outerSum is no longer increasing
            while (innerArg > eps(innerSum)) and (outerArg * (1 - innerSum) > realmin):
                m = m + 1
                k = k + 1
                innerArg = innerArg * NX / k
                innerSum = innerSum + innerArg
                outerArg = outerArg * Y / m
                outerSum = outerSum + outerArg * (1 - innerSum)

            # Form 2 also has a summation over m separate from the double sum term
            num_million = int(floor(N / 1e6))
            for m in range(0, num_million):
                outerSum = outerSum + matlabsum(expA(Y, np.array(list(range((m - 1) * int(1e6), m * int(1e6) - 1)))))

            outerSum = outerSum + sum(expA(Y, np.array(list(range(num_million * int(1e6), int(N))))))
            
            Q = outerSum
            
        else: # case 4
            # Form 4 corresponds to Eq. 11
            # Adjust start of inner summation to avoid underflow
            mmin = findStartOfSummation(Y, 0)
            
            # Adjust start of outer summation to avoid underflow
            kmin = findStartOfSummation(NX, matlabmax(mmin - N + 1, 0))
            
            # Initialize k and m to adjusted values
            k = kmin
            m = N - 1 + k
            
            # Initialize the values of innerSum and outerSum
            #
            # If value of m based on N and k is greater than the minimum value
            # computed above, pre-compute start of inner summation.  If not,
            # then we can ignore the first m-1 terms in the inner summation.
            if mmin < m:
                innerSum = 0
                for idx in range(mmin, m):
                    innerArg = expA(Y, idx)
                    if innerArg > realmin:
                        innerSum = innerSum + innerArg
            else:
                innerArg = expA(Y, m)
                innerSum = innerArg

            # With k = kmin, ignore the first k-1 terms in the outer summation
            outerArg = expA(NX, k)
            outerSum = outerArg * (1 - innerSum)
            
            # Sum until the outerSum is no longer increasing
            while (innerArg > eps(innerSum)) and (outerArg * (1 - innerSum) > realmin):
                m = m + 1
                k = k + 1
                innerArg = innerArg * Y / m
                innerSum = innerSum + innerArg
                outerArg = outerArg * NX / k
                outerSum = outerSum + outerArg * (1 - innerSum)
            
            Q = 1 - outerSum

    return Q

def findStartOfSummation(constTerm, minVal):
    epsM = 1e-40
    
    G = -4 * log(4 * epsM * (1 - epsM)) / 5 # Eq. 41

    '''
    %the parameter constTerm should be large enough for Shnidman's assumptions
    %to hold. If the parameter is not large enough, the ensuing calculations
    %are meaningless (we get a square root with negative input in Eq. 40).
    %Shnidman only offers a solution to the underflow region problem for large
    %constTerm because: 
    %1. This is effectively the case of interest 
    %2. If constTerm is not large, the solution can't be obtained in a expedient
    % manner, and is thus inefficient. 
    % Here, we assume that if (2*constTerm - G)>0, the assumptions hold and the
    % underflow region is investigated. Else, this step is skipped.
    '''
    if (expA(constTerm, minVal) > realmin) or ((2 * constTerm - G) < 0):
        # If there is no underflow in the minVal term, or if simplifying
        # assumptions do not hold, start summation at minVal
        startIdx = minVal
    else:
        # Underflow is detected in the minVal term, so compute new starting index
        startIdx = floor(constTerm + 1 / 2 - sqrt(G * (2 * constTerm - G))) # Eq. 40

    return startIdx

def expA(y, n):
    '''
     Evaluates terms of the form exp(A) = (e^-y)*(y^n)/(n!).  For large values of y
     we use a modified expression for exp(A) from a follow-up paper by Shnidman:
    
     [3] D. A. Shnidman, "Note on 'The Calculation of the Probability of
     Detection and the Generalized Marcum Q-Function'", IEEE Transactions
     on Information Theory, vol. 37, no. 4, p. 1233, July 1991.
    '''
    
    if y > 0:
        if y > 1e4:
            # Note:  There is a typo in the equation for "A" in [3], affecting the 
            # first term: (z+0.5) should read (z-0.5) 
            # and the term in the denominator of the first term in the brackets:
            # 1+1/(2z) should read 1-1/(2z).
            z = n + 1
            x = exp((z - 0.5) * ((1 - y / z) / (1 - 1 / (2 * z)) + log(y / z)) - 0.5 * log(2 * pi * y) - J(z))
        else:
            x = exp(-y + n * log(y) - gammaln(n + 1))
    else:
        # In this case y equals zero, so the answer is zero for all cases except
        # when n equals zero, in which case the answer is one
        if n == 0:
            x = 1
        else:
            x = 0

    return x

def J(z):
    # A continued fraction approximation to the Binet function given in [1]
    
    x = 1 / (12 * z + 2 / (5 * z + 53 / (42 * z + 1170 / (53 * z + 53 / z)))) # Eq. B5

    return x

def rocpfa(Pfa, SignalType = "NonfluctuatingCoherent", MinSNR = 0, MaxSNR = 20, NumPoints = 101, NumPulses = 1):
    '''
    %rocpfa Receiver operating characteristic curves on varying Pfa
    %   [Pd, SNR] = rocpfa(Pfa) returns the calculated receiver operating
    %   characteristic (ROC) curve, i.e., probability of detection (Pd) vs.
    %   single pulse signal to noise ratio (SNR) (in dB), for a given
    %   probability of false alarm (Pfa). Pd and SNR are returned in columns. 
    %
    %   If a vector of Pfa values are specified, Pd will be a matrix with each
    %   column corresponding to a given Pfa value. SNR is always a column
    %   vector and is shared by all Pfa values.
    %
    %   [Pd, SNR] = rocpfa(...,'SignalType',TYPE) specifies the type of the
    %   received signal when calculating the ROC curve. TYPE can be one of the
    %   following: 'Real' | 'NonfluctuatingCoherent' | 
    %   'NonfluctuatingNoncoherent' | 'Swerling1' | 'Swerling2' | 'Swerling3' |
    %   'Swerling4'. The default TYPE is 'NonfluctuatingCoherent'. The noise is
    %   assumed to be white Gaussian.
    %
    %   [Pd, SNR] = rocpfa(...,'MinSNR',MINSNR) specifies the minimum SNR (in
    %   dB) included in the ROC calculation. The default value of MINSNR is 0.
    %
    %   [Pd, SNR] = rocpfa(...,'MaxSNR',MAXSNR) specifies the maximum SNR (in
    %   dB) included in the ROC calculation. The default value of MAXSNR is 20.
    %
    %   [Pd, SNR] = rocpfa(...,'NumPoints',NUMPTS) specifies the number of
    %   points, i.e., steps, used to calculate the ROC curves. The default
    %   value of NUMPTS is 101.
    %
    %   [Pd, SNR] = rocpfa(...,'NumPulses',NUMINT) integrates NUMINT pulses in
    %   calculating the ROC curves. The default value of NUMINT is 1, i.e., no
    %   pulse integration.
    %
    %   rocpfa(...) plots the receiver operating characteristic curves.
    %
    %   % Example:
    %   %   Plot the ROC curve for the nonfluctuating-coherent signal type
    %   %   using a set of given Pfa values.
    %
    %   p = [1e-10 1e-8 1e-6 1e-4]; % Pfa 
    %   rocpfa(p,'SignalType','NonfluctuatingCoherent');
    %
    %   See also phased, rocsnr, npwgnthresh.
    
    %   Copyright 2010-2018 The MathWorks, Inc.
    
    %   Reference
    %   [1] Mark Richards, Fundamentals of Radar Signal Processing
    '''
    
    SNR = linspace(MinSNR, MaxSNR, NumPoints)
    d = db2pow(SNR)
    
    Pd = rocpdcalc(Pfa, d, NumPulses, SignalType)
    
    return Pd, SNR
    
def rocpdcalc(PFAs, SNRs, N, SignalType):
    '''
    % This function is for internal use and may be modified or removed in a future release.
    
    %ROCPDCALC Calculate Pd value for ROC curves
    %   Pd = ROCPDCALC(Pfa,SNR,N,TYPE) returns the probability of detection
    %   Pd corresponding the given probability of false alarm Pfa, signal to
    %   noise ratio SNR and number of pulses N used in pulse integration.
    %   TYPE specifies the signal type and can be one of the following values:
    %   [ 'Real' | 'NonfluctuatingNoncoherent' | {'NonfluctuatingCoherent'} | 
    %   'Swerling1' | 'Swerling2' | 'Swerling3' | 'Swerling4' ].
    %
    %   Both Pfa and SNR can be vectors and Pd's dimension is given by
    %   [length(Pfa) length(SNR)];
    
    %   Copyright 2008-2018 The MathWorks, Inc.
    
    %   Reference
    %   [1] Mark Richards, Fundamentals of Radar Signal Processing
    '''
    if not isinstance(PFAs, np.ndarray):
        PFAs = numpyify(PFAs)

    Yb = real(gammaincinv(1 - PFAs, N))

    match SignalType.lower():
        case "real":
            Pd = localRealROC(PFAs, SNRs, N)
        case "nonfluctuatingcoherent":
            Pd = localNonfluctuatingCoherentROC(PFAs, SNRs, N)
        case "nonfluctuatingnoncoherent":
            snr_pr = db2pow(SNRs)
            Pd = pdNonfluctuatingNoncoherent(snr_pr, N, Yb)
        case "swerling1":
            # Pd = Pd_models_main._swerling1_calculation(SNRs, np.array([N]), Yb)
            snr_pr = db2pow(SNRs)
            Pd = pdSwerling1(snr_pr, N, Yb)
        case "swerling2":
            Pd = localSwerling2ROC(PFAs, SNRs, N)
        case "swerling3":
            # calc_PD = Pd_models_main._swerling3_calculation(snr_pr, np.array([N]), Yb)
            snr_pr = db2pow(SNRs)
            Pd = pdSwerling3(snr_pr, N, Yb)
        case "swerling4":
            # calc_PD = Pd_models_main._swerling4_calculation(snr_pr, np.array([N]), Yb)
            snr_pr = db2pow(SNRs)
            Pd = pdSwerling4(snr_pr, N, Yb)
        case _:
            raise ValueError(
                f"That signal type ({SignalType}) is not recognized. Please correct it. The signal type and can be one of the following values: 'Real', 'NonfluctuatingNoncoherent', 'NonfluctuatingCoherent', 'Swerling1', 'Swerling2', 'Swerling3', or 'Swerling4', and capitalization does not matter."
            )
    
    return Pd
    
def localNonfluctuatingCoherentROC(Pfa, SNR, N):
    SNRlen = numel(SNR)
    Pd = zeros((numel(Pfa), SNRlen)) # preallocate
    
    erfcinv2PFA = erfcinv(2 * Pfa)
    for k in range(0, SNRlen):
        Pd[:, k] = 0.5 * erfc(erfcinv2PFA - (sqrt(SNR[k] * N)))
    
    return Pd
    
def localRealROC(Pfa, SNR, N):
    SNRlen = numel(SNR)
    Pd = zeros((numel(Pfa), SNRlen)) # preallocate
    
    erfcinv2PFA = erfcinv(2 * Pfa)
    SQRT2 = sqrt(2)
    for k in range(0, SNRlen):
        Pd[:, k] = 0.5 * erfc(erfcinv2PFA - (sqrt(SNR[k] * N) / SQRT2))
    
    return Pd
    
def localSwerling2ROC(Pfa, SNR, N):
    SNRlen = numel(SNR)
    Pd = zeros((numel(Pfa), SNRlen))  # preallocate
    T = real(gammaincinv(1 - Pfa, N))
    
    for k in range(0, SNRlen):
        Pd[:, k] = 1 - real(gammainc(T / (1 + SNR[k]), N))
    
    return Pd

def pdSwerling1(snr_pr, N, Yb):
    snr_db = pow2db(snr_pr)
    N = int(np.round(N))

    SNRlen = numel(snr_db)
    Pd = zeros((numel(Yb), SNRlen))  # preallocate

    if SNRlen > 1:
        for k in range(0, SNRlen):
            a = 1 / (1 + N * snr_db[k])
            
            if N == 1:
                Pd[:, k] = exp(-Yb * a)
            else:
                b = N * snr_db[k] * a
                A = repmat(exp((N - 1) * log(Yb) - Yb - gammaln(N)), 1, numel(snr_db[k]))
                Pd[:, k] = 1 - real(gammainc(Yb, N - 1)) + real(gammainc(Yb * b, N - 1)) * A
    else:
        a = 1 / (1 + N * snr_db)
            
        if N == 1:
            Pd = exp(-Yb * a)
        else:
            b = N * snr_db * a
            A = repmat(exp((N - 1) * log(Yb) - Yb - gammaln(N)), 1, numel(snr_db))
            Pd = 1 - real(gammainc(Yb, N - 1)) + real(gammainc(Yb * b, N - 1)) * A
    
    Pd[Pd > 1] = 1
    
    return Pd

def pdSwerling3(snr_pr, N, Yb):
    snr_db = pow2db(snr_pr)
    N = int(np.round(N))

    SNRlen = numel(snr_db)
    
    c = 1 / (1 + N * snr_db / 2)
    
    K0 = 1 - (N - 2) * c / (1 - c) + Yb * c
    
    if N < 3:
        Pd = K0 * exp(-Yb * c - (N - 2) * log(1 - c))
    else:
        A = repmat((N - 1) * log(Yb) - Yb, 1, SNRlen)
        Pd = exp(A + log(1 - c) - gammaln(N)) * real(gammainc(Yb * (1 - c), N - 1))
        Pd = exp(A + log(c) - gammaln(N - 1)) + 1 - real(gammainc(Yb, N - 1)) + K0 * Pd

    if SNRlen > 1:
        Pd[Pd > 1] = 1
    else:
        if Pd > 1:
            Pd = 1

    return Pd

def pdSwerling4(snr_pr, N, Yb):
    snr_db = pow2db(snr_pr)
    N = int(np.round(N))

    SNRlen = numel(snr_db)
    Pd = zeros((numel(Yb), SNRlen))

    if SNRlen > 1:
        for k in range(0, SNRlen):
            tmp1 = Yb / (1 + snr_db[k] / 2)
            pdtmp = zeros((numel(Yb), 1))
            for gammaI in range(0, N):
                c = gammaI * log(snr_db[k] / 2) + sum(log(N - gammaI + np.arange(1, N))) - sum(log(np.arange(1, gammaI)))
                pdtmp = pdtmp + exp(c - N * log(1 + snr_db[k] / 2)) * real(gammainc(tmp1, N + gammaI))
    
            Pd[:, k] = 1 - pdtmp
    else:
        tmp1 = Yb / (1 + snr_db / 2)
        pdtmp = zeros((numel(Yb), 1))
        for gammaI in range(0, N + 1):
            c = gammaI * log(snr_db / 2) + sum(log(N - gammaI + np.arange(1, N))) - sum(log(np.arange(1, gammaI)))
            pdtmp = pdtmp + exp(c - N * log(1 + snr_db / 2)) * real(gammainc(tmp1, N + gammaI))

        Pd = 1 - pdtmp
    
    Pd[Pd > 1] = 1

    return unlevel(Pd)

def pdNonfluctuatingNoncoherent(snr_pr, N, Yb):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %   Copyright 2020 The MathWorks, Inc.
    '''

    snr_db = pow2db(snr_pr)

    SNRlen = numel(snr_db)
    Pd = zeros((numel(Yb), SNRlen)) # preallocate
    # T = real(gammaincinv(1 - Pfa, N))
    # T = T(:)
    N = np.round(N).astype(int)
    r = np.arange(1, N)
    rr = repmat(r, numel(Yb), 1)

    if SNRlen > 1:
        for k in range(0, SNRlen):
            NX = N * snr_db[k]
            Pd[:, k] = marcumq(sqrt(2 * NX), sqrt(2 * Yb))
            SQRTNXT2 = 2 * sqrt(NX * Yb)
            if (snr_db[k] != 0) and (snr_db[k] != inf) and (numel(r) > 0):
                A = bsxfun("minus", log(Yb / NX) * r / 2, ((Yb + NX) - SQRTNXT2))
                Pdtemp = exp(A) * real(besseli(rr, repmat(SQRTNXT2, 1, N - 1), 1))
                Pdtemp[~isfinite(Pdtemp)] = 0
                Pd[:, k] = Pd[:, k] + matlabsum(Pdtemp, 2)
    else:
        NX = N * snr_db
        Pd = marcumq(sqrt(2 * NX), sqrt(2 * Yb))
        SQRTNXT2 = 2 * sqrt(NX * Yb)
        if (snr_db != 0) and (snr_db != inf) and (numel(r) > 0):
            A = bsxfun("minus", log(Yb / NX) * r / 2, ((Yb + NX) - SQRTNXT2))
            Pdtemp = exp(A) * real(besseli(rr, repmat(SQRTNXT2, 1, N - 1), 1))
            Pdtemp[~isfinite(Pdtemp)] = 0
            Pd = Pd + matlabsum(Pdtemp, 2)

    Pd[Pd > 1] = 1

    return unlevel(Pd)

def radarmetricplot(
    r,
    metric,
    objective = None,
    threshold = None,
    MaxRangeRequirement = None,
    ShowStoplight = False,
    RadarName = None,
    MetricName = "Metric",
    RequirementName = "Requirement",
    RangeUnit = 'm',
    MetricUnit = '',
    Parent = None,
):
    # Create new axis object if not specified
    if Parent is None:
        fig, Parent = plt.subplots(figsize=(8, 5))

    # Horizontal plot limits
    RangeLims = [r[0], r[-1]]
    
    # Vertical plot limits
    ydataMinvals = []
    ydataMaxvals = []
    
    ydataMinvals.append(matlabmin(np.asarray(metric)[~isinf(metric)], [], 'all'))
    ydataMaxvals.append(matlabmax(np.asarray(metric)[~isinf(metric)], [], 'all'))
    
    if objective is not None:
        ydataMinvals.append(matlabmin(np.asarray(objective)[~isinf(objective)], [], 'all'))
        ydataMaxvals.append(matlabmax(np.asarray(objective)[~isinf(objective)], [], 'all'))
    if threshold is not None:
        ydataMinvals.append(matlabmin(np.asarray(threshold)[~isinf(threshold)], [], 'all'))
        ydataMaxvals.append(matlabmax(np.asarray(threshold)[~isinf(threshold)], [], 'all'))
    
    ydataMin, _ = matlabmin(ydataMinvals)
    ydataMax, _ = matlabmax(ydataMaxvals)
    ydataSpan = ydataMax - ydataMin
    YLim = [ydataMin - ydataSpan * 0.1, ydataMax + ydataSpan * 0.1]

    # Create a label for each radar
    Radar_labels = []
    if RadarName is None: # if RadarName is not specified
        if np.ndim(r) > 1 and np.ndim(metric) > 1: # if more than one radar is described
            for Radar_idx in range(np.shape(r)[1]): # Loop through radars
                Radar_labels.append(f"Radar{Radar_idx + 1}")
                Parent.plot(r[Radar_idx], metric[Radar_idx], label = Radar_labels[Radar_idx], linewidth = 2)
        elif np.ndim(r) == 1 and np.ndim(metric) == 1: # if more than one radar is described
            Radar_labels.append("Radar1")
            Parent.plot(r, metric, color = 'b', label = Radar_labels[0], linewidth = 2)
        else:
            raise ValueError("r and/or metric dimensions are inappropriate")
    else:
        Radar_labels = RadarName
                
        if np.ndim(r) > 1 and np.ndim(metric) > 1: # if more than one radar is described
            assert np.shape(r)[1] == numel(RadarName)
            
            for Radar_idx in range(np.shape(r)[1]): # Loop through radars
                Radar_labels.append(f"Radar{Radar_idx + 1}")
                Parent.plot(r[Radar_idx], metric[Radar_idx], label = Radar_labels[Radar_idx], linewidth = 2)
        elif np.ndim(r) == 1 and np.ndim(metric) == 1: # if more than one radar is described
            Parent.plot(r, metric, color = 'b', label = Radar_labels, linewidth = 2)
        else:
            raise ValueError("r and/or metric dimensions are inappropriate")

    if (MaxRangeRequirement is not None) and (MaxRangeRequirement >= RangeLims[0]) and (MaxRangeRequirement <= RangeLims[1]):
        Parent.axvline(x = MaxRangeRequirement, color = 'k', label = 'Max Range', linestyle = '--', linewidth = 2)

        # Transform data coordinates to display (pixel) coordinates and then to axes fraction coordinates
        # Parent.transAxes.inverted() creates a transform that goes from display (pixel)
        # back to the axes fraction coordinates (where (0,0) is bottom-left, (1,1) is top-right)
        display_coords = Parent.transData.transform((MaxRangeRequirement, 0))
        MaxRangeRequirement_display_coords = Parent.transAxes.inverted().transform(display_coords)
    else:
        MaxRangeRequirement_display_coords = 0

    if objective is not None:
        Parent.axhline(y = objective, color = 'b', label = RequirementName, linestyle = '--', linewidth = 2)
        if ShowStoplight:
            # Highlight a horizontal 'pass' region
            # Y is in data units, x is in display units
            Parent.axhspan(objective, YLim[1], facecolor = 'limegreen', alpha = .5)
            if threshold is None:
                Parent.axhspan(YLim[0], objective, facecolor = 'red', alpha = .5)

            if (MaxRangeRequirement is not None) and (MaxRangeRequirement >= RangeLims[0]) and (MaxRangeRequirement <= RangeLims[1]):    
                # Highlight a horizontal 'grey' region
                # Y is in data units, x is in display units
                Parent.axhspan(YLim[0], objective, xmin = MaxRangeRequirement_display_coords[0], xmax = 1, facecolor = 'grey', alpha = 1)

    if threshold is not None:
        Parent.axhline(y = threshold, color = 'b', label = 'Threshold', linestyle = '-.', linewidth = 2)
        if ShowStoplight and (objective is not None):
            # Highlight a horizontal 'warn' region
            # Y is in data units, x is in display units
            Parent.axhspan(threshold, objective, facecolor = 'yellow', alpha = .5)

            # Highlight a horizontal 'fail' region
            # Y is in data units, x is in display units
            Parent.axhspan(YLim[0], threshold, facecolor = 'red', alpha = .5)

            if (MaxRangeRequirement is not None) and (MaxRangeRequirement >= RangeLims[0]) and (MaxRangeRequirement <= RangeLims[1]):    
                # Highlight a horizontal 'grey' region
                # Y is in data units, x is in display units
                Parent.axhspan(YLim[0], threshold, xmin = MaxRangeRequirement_display_coords[0], xmax = 1, facecolor = 'grey', alpha = 1)
            
        elif ShowStoplight: # objective is None
            # Highlight a horizontal 'warn' region
            # Y is in data units, x is in display units
            Parent.axhspan(threshold, YLim[1], facecolor = 'yellow', alpha = .5)

            # Highlight a horizontal 'fail' region
            # Y is in data units, x is in display units
            Parent.axhspan(YLim[0], threshold, facecolor = 'red', alpha = .5)

            if (MaxRangeRequirement is not None) and (MaxRangeRequirement >= RangeLims[0]) and (MaxRangeRequirement <= RangeLims[1]):    
                # Highlight a horizontal 'grey' region
                # Y is in data units, x is in display units
                Parent.axhspan(YLim[0], threshold, xmin = MaxRangeRequirement_display_coords[0], xmax = 1, facecolor = 'grey', alpha = 1)
    
    Parent.set_xlabel(f"Target Range ({RangeUnit})")
    Parent.set_ylabel(f"{MetricName} ({MetricUnit})")
    Parent.set_title(f'{MetricName} vs. Range')
    Parent.legend()
    Parent.grid(True)
    Parent.set_xlim(RangeLims)
    Parent.set_ylim(YLim)

    return Parent

def convert2meters(fromUnit):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %convert2meters Conversion factor from input unit to meters. 
    %   FAC = convert2meters(FROM) converts FROM unit to meters. Accepted input
    %   units are 'nmi' | 'mi' | 'km' | 'ft' | 'm' | 'kft'.
    % 
    %   See the calling function for additional information.
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also radarvcd, unitsratio.
    
    %   Copyright 2021 The MathWorks, Inc.
    '''
    
    conversionTableFromUnit = ['nmi', 'nm', 'mi', 'km', 'ft', 'kft']
    conversionTableFrom = [1852, 1852, 5280 * 0.3048, 1e3, 0.3048, 1e3 * 0.3048]
    numConversions = numel(conversionTableFrom)
    fac = 1
    notFound = True
    ii = 0

    while notFound and (ii <= numConversions):
        if conversionTableFromUnit[ii] == fromUnit:
            fac = conversionTableFrom[ii]
            notFound = False
        ii = ii + 1

    return fac

def radarvcd(
    freq, # Radar frequency
    rfs, # Free-space range
    anht, # Radar antenna height
    rngunit = "km", # RangeUnit
    htunit = "m", # HeightUnit
    epsc = None, # SurfaceRelativePermittivity
    pol = "H", # Polarization
    anpat = None,  # AntennaPattern
    patel = None,  # PatternAngles
    tilt = 0, # TiltAngle
    htsd = 0, # SurfaceHeightStandardDeviation
    beta0 = 0, # SurfaceSlope
    vegType = "None", # VegetationType
    elbw = 10, # ElevationBeamwidth
    Re = Rearth, # EffectiveEarthRadius
    vcpelmax = 60, # MaxElevation
    vcpelmin = 0, # MinElevation
    thdel = 1, # ElevationStepSize
):
    '''
    %radarvcd Vertical coverage diagram for radar
    %   [VCP, VCPANGLES] = radarvcd(FREQ,RFS,ANHT) outputs the vertical
    %   coverage pattern VCP at elevation angles VCPANGLES (in degrees) for a
    %   radar system operating at a positive, scalar carrier frequency FREQ (in
    %   Hz). RFS is the calculated or assumed free-space range for a target or
    %   for a one-way RF system at which the field strength would have a
    %   specified value. RFS is specified as an N-length vector of positive
    %   values. ANHT is the nonnegative, scalar antenna height, assumed to be
    %   referenced to the surface. By default, the unit of RFS is km, and the
    %   unit of ANHT is m.
    %
    %   The output VCP is an MxN matrix where N corresponds to the number of
    %   ranges specified in RFS. Each column represents an individual vertical
    %   coverage contour. The default units of VCP is km. VCPANGLES is an
    %   M-length column vector. Each entry in VCPANGLES specifies the elevation
    %   angle (in degrees) at which the vertical coverage pattern is measured.
    %
    %   The vertical coverage pattern is a detection or constant-signal-level
    %   contour in the vertical plane that takes into account the interference
    %   between the direct and ground-reflected rays. Atmospheric refraction is
    %   modeled through the use of an effective Earth radius. Scattering and
    %   ducting are assumed to be negligible.
    %
    %   The vertical coverage pattern is generally considered to be valid for
    %   antenna heights that are within a few hundred feet of the surface and
    %   with targets at altitudes that are not too close to the radar horizon.
    %   
    %   [...] = radarvcd(...,'RangeUnit',RNGUNIT) specifies the units of range
    %   RFS as one of 'nmi' | 'mi' | 'km' | 'ft' | 'm' | 'kft'. The default
    %   value is 'km'.
    % 
    %   [...] = radarvcd(...,'HeightUnit',HTUNIT) specifies the units of height
    %   ANHT as one of 'nmi' | 'mi' | 'km' | 'ft' | 'm' | 'kft'. The default
    %   value is 'm'.
    % 
    %   [...] = radarvcd(...,'Polarization',pol) specifies the polarization of
    %   the transmitted wave. The parameter, pol, can be one of 'H' | 'V',
    %   where the default value is 'H'. 'H' indicates horizontal polarization
    %   and 'V' indicates vertical polarization.
    % 
    %   [...] = radarvcd(...,'SurfaceRelativePermittivity',EPSC) specifies the
    %   scalar, complex permittivity (dielectric constant) of the reflecting
    %   surface. The default value of EPSC depends on the value of FREQ. The
    %   function uses a sea water model that is valid up to 10 GHz.
    % 
    %   [...] = radarvcd(...,'SurfaceHeightStandardDeviation',HGTSD) specifies
    %   the scalar standard deviation of the surface height in meters. The
    %   default value of HGTSD is 0, indicating a smooth surface. The unit of
    %   height is specified by HTUNIT.
    %
    %   [...] = radarvcd(...,'SurfaceSlope',BETA0) specifies the nonnegative
    %   scalar surface slope in degrees. This value is expected to be 1.4 times
    %   the RMS surface slope. Given the condition that
    %        2*GRAZ/BETA0 < 1,
    %   where GRAZ is the grazing angle of the geometry specified in degrees,
    %   the effective surface height standard deviation in meters is calculated
    %   as
    %        Effective HGTSD = HGTSD*(2*GRAZ/BETA0)^(0.2),
    %   which better accounts for shadowing. Otherwise, the effective height
    %   standard deviation is equal to HGTSD. BETA0 defaults to 0 degrees,
    %   indicating a smooth surface.
    %   
    %   [...] = radarvcd(...,'VegetationType',VEGTYPE) specifies the vegetation
    %   type of the surface as one of either 'Trees' | 'Weeds' | 'Brush' |
    %   'Grass' | 'None'. In the case of VEGTYPE set to 'Trees' | 'Weeds' |
    %   'Brush', there is an assumption of dense vegetation. The case of
    %   VEGTYPE set to 'Grass' assumes thin grass. Use this argument when using
    %   the function on surfaces different from the sea. Defaults to 'None'.
    % 
    %   [...] = radarvcd(...,'ElevationBeamwidth',ELBW) specifies the positive
    %   scalar half-power elevation beamwidth in degrees. The elevation
    %   beamwidth is used in the calculation of a sinc antenna pattern. The
    %   default antenna pattern is symmetrical with respect to the beam maximum
    %   and is of the form sin(u)/u. The parameter u is given by
    %        u = k*sind(theta), 
    %   where theta is the elevation angle in degrees and k is given by
    %        k = 1.39157/sind(ELBW/2). 
    %   ELBW defaults to 10 degrees.
    %
    %   [...] = radarvcd(...,'AntennaPattern',PAT,'PatternAngles', PATEL)
    %   specifies the antenna elevation pattern's normalized voltage response
    %   in linear units (Volts) and corresponding elevation angles in degrees.
    %   This is an alternative to specifying the elevation beamwidth. Both PAT
    %   and PATEL must be vectors of the same size. PATEL must be between -90
    %   and 90 degrees. In general, to properly compute the coverage, the
    %   pattern should be specified from -90 to 90 degrees. If both an antenna
    %   pattern and an elevation beamwidth are provided, the function uses the
    %   antenna pattern and ignores the elevation beamwidth value. Defaults to
    %   a sinc antenna pattern.
    % 
    %   [...] = radarvcd(...,'TiltAngle',TILTANG) specifies the scalar tilt
    %   angle in degrees of the antenna with respect to the surface. TILTANG
    %   defaults to 0 degrees.
    %
    %   [...] = radarvcd(...,'EffectiveEarthRadius',RE) specifies the
    %   effective Earth radius as a positive scalar in meters. The effective
    %   Earth radius is an approximation used for modeling refraction effects
    %   in the troposphere. The default calculates the effective Earth radius
    %   using a refractivity gradient of -39e-9, which results in approximately
    %   4/3 of the real Earth radius.
    % 
    %   [...] = radarvcd(...,'MaxElevation', MAXEL) specifies the maximum
    %   elevation angle, in degrees, for which the vertical coverage pattern
    %   will be calculated. The default value of MAXEL is 60.
    %
    %   [...] = radarvcd(...,'MinElevation', MINEL) specifies the nonnegative
    %   minimum elevation angle, in degrees, for which the vertical coverage
    %   pattern will be calculated. The default value of MINEL is 0.
    %
    %   [...] = radarvcd(...,'ElevationStepSize', ELSTEP) specifies the
    %   elevation angle increment value as a positive scalar, in degrees, for
    %   creating the elevation vector, which is specified as
    %   MINEL:ELSTEP:MAXEL. The default value of ELSTEP is calculated as
    %       ELSTEP = 885.6/(pi*freqMHz*anhtFt),
    %   where freqMHz is the frequency in MHz and anhtFt is the antenna height
    %   in feet.
    %
    %   radarvcd(...) plots the radar coverage diagram.
    %
    %   % Example:
    %   %   Plot the radar vertical coverage pattern assuming the antenna has a
    %   %   sinc pattern. The frequency is 100 MHz, the antenna height is 20
    %   %   feet, and the range is 100 nautical miles. Assume the surface is
    %   %   smooth, the antenna is not tilted, and the transmitted polarization
    %   %   is horizontal.
    %
    %   patAng     = linspace(-90,90,361)';
    %   pat_u      = 1.39157/sind(90/2)*sind(patAng);
    %   pat        = sinc(pat_u/pi);
    % 
    %   freq      = 100e6; % Frequency (Hz)
    %   anht      = 20;    % Antenna height (ft)
    %   rfs       = 100;   % Range (nmi) 
    %   tiltAng   = 0;     % Tilt angle (deg) 
    %   radarvcd(freq,rfs,anht,'RangeUnit','nmi','HeightUnit','ft', ...
    %       'AntennaPattern',pat,'PatternAngles',patAng,'TiltAngle',tiltAng);
    %
    %   See also radar, blakechart, effearthradius.
    
    %   Copyright 2012-2021 The MathWorks, Inc.
    
    %   References
    %     [1] Blake, L.V. "Machine Plotting of Radar Vertical-Plane Coverage
    %         Diagrams." Naval Research Laboratory, NRL Report 7098, 1970.
    %     [2] Barton, David K. Radar Equations for Modern Radar. Norwood, MA:
    %         Artech House, 2013.
    
    %#ok<*EMCLS>
    %#ok<*EMCA>
    %#codegen
    '''
    
    # Calculate elevation step size
    if thdel == 1:
        # Code snippet from NRL Report 7098, p38 
        if anht < sqrt(eps()):
            thdel = pi / 1800 # Guards against case when anht is 0 and thdel becomes inf
        else:
            m2ft = 1 / 0.3048
            thdel = 4.92 / m2ft / (freq * 1e-6 * anht) # Elevation step size (rad)
    else:
        thdel = deg2rad(thdel)

    if (anpat is None) or (patel is None):
        anpat, patel = sincpattern(elbw)

    if epsc is None:
        epsc = seaComplexPermittivity(freq)
    
    # Sort free space range and make sure it is a row vector 
    rfsC = sort(rfs)

    # Convert angles to radians
    vcpelmax = deg2rad(vcpelmax)
    vcpelmin = deg2rad(vcpelmin)
    
    # Calculate direct ray elevation angles
    thetad = np.arange(vcpelmin, vcpelmax, thdel) # Direct ray elevation angles

    if htunit != "m":
        htfac = convert2meters(htunit)
    else:
        htfac = 1
        
    # Convert anht and htsd to meters
    anht = anht * htfac
    htsd = htsd * htfac

    # assert anht < htsd

    # Calculate propagation factor 
    pf, pd, wv = propfactorinterf(freq, anht, thetad, epsc, pol, anpat, patel, tilt, htsd, beta0, vegType, Re)

    # Limit theta
    angIdx = find(pd >= wv / 4, 1)
    startIdx = 1
    if ~isempty(angIdx):
        startIdx = angIdx[0]
    
    # Calculate vertical coverage pattern, same unit as input range RFS
    # VCP is a matrix with each column representing a separate vertical
    # coverage contour
    vcp = pf[startIdx:] * rfsC
    vcp[isnan(vcp)] = 0
    vcpangles = rad2deg(thetad[startIdx:])
    
    # # Output results
    # if ~nargout % Output args:
    #     if ~isempty(coder.target):
    #         coder.internal.assert(false,'phased:rocsnr:invalidCodegenOutput','radarvcd')

    #     blakechart(vcp, vcpangles,'RangeUnit', rngunit, 'HeightUnit', htunit, 'ScalePower', 1)

    return vcp, vcpangles

def blakechart(
    vcp, # Vertical coverage pattern
    vcpangles, # Vertical coverage pattern angles
    rmax = None, # Maximum range of plot
    hmax = None, # Maximum height of plot
    rngunit = 'km', # RangeUnit
    htunit = 'km', # HeightUnit
    spow = 0.25, # ScalePower
    Ns = 313, # SurfaceRefractivity
    expo = 0.143859, # RefractionExponent
    anht = 0, # AntennaHeight
    fc = 'blue', # FaceColor
    ec = 'blue', # EdgeColor
    hAxes = None, # Parent axes
):
    '''
    %blakechart Range-angle-height (Blake) chart
    %   blakechart(VCP,VCPANGLES) plots the vertical coverage pattern, VCP, for
    %   a radar system on a Blake chart, also known as a range-height-angle
    %   chart. The vertical coverage pattern is a detection or
    %   constant-signal-level contour. Normal atmospheric refraction is taken
    %   into account through the use of the CRPL exponential reference
    %   atmosphere with a refractivity of 313 N-units and a refraction exponent
    %   (decay constant) of 0.143859 1/km. Scattering and ducting are assumed
    %   to be negligible.
    %   
    %   The range in the range-height-angle chart is the propagated range and
    %   the height is relative to the origin of the ray. It is assumed that the
    %   antenna height is not at an appreciable height above ground level (<
    %   1000 ft or about 305 m).
    %
    %   VCP can be a matrix whose columns represent individual vertical
    %   coverage patterns. VCPANGLES is a column vector whose number of rows is
    %   the same as number of rows of VCP. Each entry in VCPANGLES specifies
    %   the elevation angle (in degrees) at which the vertical coverage pattern
    %   is measured.
    %
    %   blakechart(VCP,VCPANGLES,RMAX,HMAX) specifies the range limit RMAX and
    %   height limit HMAX for the Blake chart, respectively. The units of RMAX
    %   and HMAX are determined by the values of the RangeUnit and HeightUnit
    %   parameters. The default unit for both is 'km'.
    % 
    %   blakechart(...,'RangeUnit',RNGUNIT,'HeightUnit',HTUNIT) specifies the
    %   range and height units in RNGUNIT and HTUNIT as one of 'km' | 'm' |
    %   'mi' | 'nmi' | 'ft' | 'kft', where the default values for both are
    %   'km'.
    % 
    %   blakechart(...,'ScalePower',SPOW) specifies the range and height axis
    %   scale power as a scalar between 0 and 1. The default value for SPOW is
    %   1/4.
    % 
    %   blakechart(...,'SurfaceRefractivity',NS,'RefractionExponent',REXP)
    %   specifies the surface refractivity NS (in N-Units) as a non-negative
    %   scalar and the refraction exponent factor (decay constant) REXP (in
    %   1/km) as a non-negative scalar for the CRPL exponential reference
    %   atmosphere model. The default value of NS is 313 and the default value
    %   of REXP is 0.143859.
    %
    %   The atmospheric refraction model is CRPL exponential reference model
    %   given by
    %        n(h) = 1 + NS*1e-6*exp(-REXP*h),
    %   where h is the height in kilometers and n is the index of refraction.
    %
    %   blakechart(...,'AntennaHeight',ANHT) specifies the antenna height. The
    %   unit of height is specified by HTUNIT. When an antenna height is
    %   provided, the height in the Blake chart is the height above ground
    %   level. Otherwise, the height in the Blake chart is relative to the
    %   origin of the ray, and it is assumed that the antenna height is not at
    %   an appreciable height above ground level (< 1000 ft or about 305 m).
    %   Default antenna height is 0.
    %
    %   blakechart(...,'FaceColor',FC) specifies the face color of the vertical
    %   coverage pattern patch. FC can be specified by the color name (e.g.,
    %   'red'), the short name ('r'), a hexadecimal color code ('#FF0000'), or
    %   an RGB triplet ([1 0 0]). If specifying more than one color, the number
    %   of colors must match the number of columns in the VCP input. FC can
    %   also be set to 'none' if no patch fill is desired. FC defaults to the
    %   default MATLAB color order. See MATLAB ColorSpec for additional
    %   information.
    %
    %   blakechart(...,'EdgeColor',EC) specifies the edge color of the vertical
    %   coverage pattern patch. EC can be specified by the color name (e.g.,
    %   'red'), the short name ('r'), a hexadecimal color code ('#FF0000'), or
    %   an RGB triplet ([1 0 0]). If specifying more than one color, the number
    %   of colors must match the number of columns in the VCP input. EC can
    %   also be set to 'none' if no edge color is desired. EC defaults to the
    %   default MATLAB color order. See MATLAB ColorSpec for additional
    %   information.
    %
    %   blakechart(...,'Parent',HAX) specifies the plot axes, HAX. The default
    %   axes are in the current figure.
    % 
    %   % Example:
    %   %   Plot the radar vertical coverage pattern assuming the antenna has a
    %   %   sinc pattern. The frequency is 100 MHz, the antenna height is 20
    %   %   feet, and the range is 100 nautical miles. Assume the surface is 
    %   %   smooth, the antenna is not tilted, and the transmitted polarization
    %   %   is horizontal.
    %
    %   patAng    = linspace(-90,90,361)';
    %   pat_u     = 1.39157/sind(90/2)*sind(patAng);
    %   pat       = sinc(pat_u/pi);
    % 
    %   freq      = 100e6;
    %   anht      = 20;
    %   rfs       = 100;
    %   tiltAng   = 0;
    %   [vcdNmi, vcdAng] = radarvcd(freq,rfs,anht, ...
    %       'RangeUnit','nmi','HeightUnit','ft', ...
    %       'AntennaPattern',pat,'PatternAngles',patAng,'TiltAngle',tiltAng);
    %
    %   blakechart(vcdNmi, vcdAng, ...
    %       'RangeUnit','nmi','HeightUnit','ft')
    %   
    %   See also radar, radarvcd, refractionexp.
    
    %   Copyright 2012-2023 The MathWorks, Inc.
    
    %   References:
    %
    %   [1] Blake, L.V. "Radio Ray (Radar) Range-Height-Angle Charts." Naval
    %       Research Laboratory, NRL Report 6650, Jan. 22, 1968.
    %   [2] Blake, L.V. "Ray Height Computation for a Continuous Nonlinear
    %       Atmospheric Refractive-Index Profile." RADIO SCIENCE, Vol. 3 (New
    %       Series), No. 1, Jan. 1968, pp. 85-92.
    %   [3] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
    %       Atmosphere. Washington, DC: US Gov. Print. Off., 1959.
    '''
    userSetMaxLims = False
    if (rmax is not None) and (hmax is not None):
        userSetMaxLims = True
        
    rngunit, rngunitDisp, rngfac = unitnamemod(rngunit)
    htunit, htunitDisp, htfac = unitnamemod(htunit)
    addDataTips = False
    
    hVCP = rhaplot(userSetMaxLims, vcp, vcpangles, rngunit, rngunitDisp, rngfac, htunit, htunitDisp, htfac, spow, Ns, expo, anht, fc, ec, hAxes, addDataTips, rmax, hmax)

def unitsratio(to_in, from_in):
    '''
    %UNITSRATIO Unit conversion factors
    %
    %   RATIO = UNITSRATIO(TO, FROM) returns the number of TO units per one
    %   FROM unit.  For example, UNITSRATIO('cm', 'm') returns 100 because
    %   there are 100 centimeters per meter.  UNITSRATIO makes it easy to_in
    %   convert from_in one system of units to_in another.  Specifically, if X is
    %   in units FROM and Y is in units TO, then the conversion equation is
    %
    %                  Y = UNITSRATIO(TO, FROM) * X.
    %
    %   TO and FROM may be any of the length units supported by
    %   validateLengthUnit, or may be one of the following angle units:
    %
    %     Angle Unit            Valid Inputs
    %     ----------            ------------
    %     radian               'rad', 'radian', 'radians'
    %     degree               'deg', 'degree', 'degrees'
    %
    %   Examples
    %   --------
    %   % Approximate mean earth radius in meters
    %   radiusInMeters = earthRadius('meters')
    %   % Conversion factor
    %   feetPerMeter = unitsratio('feet', 'meter')
    %   % Radius in (international) feet:
    %   radiusInFeet = feetPerMeter * radiusInMeters
    %
    %   % The following prints a true statement for any valid TO, FROM pair:
    %   to_in   = 'feet';
    %   from_in = 'mile';
    %   sprintf('There are %g %s per %s.', unitsratio(to_in,from_in), to_in, from_in)
    %
    %   % The following prints a true statement for any valid TO, FROM pair:
    %   to_in   = 'degrees';
    %   from_in = 'radian';
    %   sprintf('One %s is %g %s.', from_in, unitsratio(to_in,from_in), to_in)
    %
    %   See also validateLengthUnit.
    
    % Copyright 2002-2017 The MathWorks, Inc.
    '''
    # Validate units and convert to standard names.
    degreeStrings = ['deg', 'degree', 'degrees']
    radianStrings = ['rad', 'radian', 'radians']
    
    try:
        to_in = validateLengthUnit(to_in, 'UNITSRATIO', 'TO', 1)
    except Exception as e:
        if any(strcmpi(to_in, degreeStrings)):
            to_in = 'degree'
        elif any(strcmpi(to_in, radianStrings)):
            to_in = 'radian'
        else:
            raise e
    
    try:
        from_in = validateLengthUnit(from_in, 'UNITSRATIO', 'FROM', 2)
    except:
        if any(strcmpi(from_in, degreeStrings)):
            from_in = 'degree'
        elif any(strcmpi(from_in, radianStrings)):
            from_in = 'radian'
        else:
            exception.throw()
    
    # Define a relationship graph by specifying a scaling factor for each pair
    # of adjacent units (using the standard names). Each unit on the left is
    # defined in terms of the unit on the right via the factor provided, which
    # is the value that will be returned by unitsratio(RIGHT, LEFT). For
    # example, unitsratio('meter','micron') returns 1e-6.
    graph = [
            ['micron',             'meter',   1e-6],
            ['millimeter',         'meter',   1e-3],
            ['centimeter',         'meter',   1e-2],
            ['kilometer',          'meter',   1e+3],
            ['nautical mile',      'meter',   1852],
            ['foot',               'meter',   0.3048],
            ['mile',               'foot',    5280],
            ['inch',               'foot',    1/12],
            ['yard',               'foot',    3],
            ['U.S. survey foot',   'meter',   1200/3937],
            ['U.S. survey mile',   'U.S. survey foot', 5280],
            ['Clarke''s foot',     'meter',   0.3047972654],
            ['German legal metre', 'meter',   1.0000135965],
            ['Indian foot',        'meter',   12/39.370142],
            ['degree',             'radian',  pi/180],
            ]
    graph = pd.DataFrame(graph)
    
    # Do a depth-first search of the directed graph corresponding
    # to_in the definitions array, recursively searching for a
    # path from_in FROM to_in TO.
    ratio, _ = searchgraph(to_in, from_in, graph, [])

    # A return value of NaN in RATIO indicates that no connection exists.
    if isnan(np.array(ratio)):
        raise ValueError(f"unitsratio: unable to_in convert units to_in {to_in} from_in {from_in}")

    return unlevel(ratio)
        
def searchgraph(to_in, from_in, graph, history):
    
    # Assume a dead-end unless/until a path is found from_in FROM to_in TO.
    ratio = nan
    
    # Stop here if FROM has already been checked (avoid loops in the graph).
    if np.any(strcmp(from_in, history)):
        return ratio, history
    
    # Append FROM to_in the list of nodes that have been visited.
    if numel(history) > 0:
        history[-1] = from_in
    
    # Find occurrences of FROM and TO in columns 1 and 2 of GRAPH.
    from1 = find(strcmp(from_in, graph.iloc[:, 0].values))
    from2 = find(strcmp(from_in, graph.iloc[:, 1].values))
    to1   = find(strcmp(to_in,   graph.iloc[:, 0].values))
    to2   = find(strcmp(to_in,   graph.iloc[:, 1].values))
    
    # See if there's a direct conversion from_in TO to_in FROM:
    # If there's a row with TO in column 1 and FROM in column 2, then
    # column 3 of that row contains the conversion factor
    # from_in FROM to_in TO.
    i = intersect(to1, from2)
    if numel(i) == 1:
        ratio = 1 / graph.iloc[i, 2]
        return ratio, history
    
    # See if there's a direct conversion from_in FROM to_in TO:
    # If there's a row with FROM in column 1 and TO in column 2,
    # then column 3 of that row contains the conversion factor.
    i = intersect(to2, from1)
    if numel(i) == 1:
        ratio = graph.iloc[i, 2]
        return ratio, history
    
    # Recursively search for conversion to_in TO from_in each node adjacent
    # to_in FROM.
    
    # Search from_in the adjacent nodes with a direct conversion _from_
    # FROM.  If a conversion factor (non-NaN) to_in TO is found from_in
    # one of these adjacent nodes, then multiply it by the conversion
    # factor from_in FROM to_in that neighbor (divide by the defining
    # factor in column 3 of GRAPH).
    for i in range(0, numel(from2)):
       n = from2[i]
       ratio, history = searchgraph(to_in, graph.iloc[n, 0], graph, history);
       if ~isnan(ratio):
           ratio = ratio / graph.iloc[n, 2]
           return ratio, history
    
    # Search from_in the adjacent nodes with a direct conversion _to_ FROM.
    # If a conversion factor (non-NaN) to_in TO is found from_in one of these
    # adjacent nodes, then divide it by the conversion factor from_in FROM
    # to_in that neighbor (multiply by the defining factor in column 3).
    for i in range(0, numel(from1)):
       n = from1[i]
       ratio, history = searchgraph(to_in, graph.iloc[n, 1], graph, history)
       if ~isnan(ratio):
           ratio = ratio * graph.iloc[n, 2]
           return ratio, history
    
    return ratio, history

def validateLengthUnit(unit, funcName = None, varName = None, argIndex = None):
    '''
    %validateLengthUnit Validate and standardize length unit name
    %
    %   standardName = validateLengthUnit(UNIT, FUNCNAME, VARNAME, ARGINDEX)
    %   checks that UNIT is a valid length unit and converts it to a standard
    %   unit name.  The following table lists each standard name that is
    %   supported, and the values that correspond to that name. The function is
    %   case-insensitive with respect to its input.  Spaces, periods, and
    %   apostrophes are ignored.  Plural forms are accepted in most cases, but
    %   the result (standardName) is always singular. The optional inputs
    %   FUNCNAME, VARNAME, and ARGINDEX may be included for use in error
    %   message formatting, with behavior identical to that provided by the
    %   VALIDATEATTRIBUTES inputs of the same names.
    %
    %     Standard Name            Supported Names
    %     -------------            ---------------
    %     meter                'm', 'meter(s)', 'metre(s)'
    %
    %     centimeter           'cm', 'centimeter(s)', 'centimetre(s)'
    %
    %     millimeter           'mm', 'millimeter(s)', 'millimetre(s)'
    %
    %     micron               'micron(s)'
    %
    %     kilometer            'km', 'kilometer(s)', 'kilometre(s)'
    %
    %     nautical mile        'nm', 'naut mi', 'nautical mile(s)'
    %
    %     foot                 'ft',   'international ft',
    %                          'foot', 'international foot',
    %                          'feet', 'international feet'
    %
    %     inch                 'in', 'inch', 'inches'
    %
    %     yard                 'yd', 'yds', 'yard(s)'
    %
    %     mile                 'mi', 'mile(s)', 'international mile(s)'
    %
    %     U.S. survey foot     'sf',
    %                          'survey ft',   'US survey ft', 'U.S. survey ft',
    %                          'survey foot', 'US survey foot', 'U.S. survey foot',
    %                          'survey feet', 'US survey feet', 'U.S. survey feet',
    %
    %     U.S. survey mile     'sm', 'survey mile(s)', 'statute mile(s)',
    %     (statute mile)       'US survey mile(s)', 'U.S. survey mile(s)'
    %
    %     Clarke's foot        'Clarke''s foot', 'Clarkes foot'
    %
    %     German legal metre   'German legal metre', 'German legal meter'
    %
    %     Indian foot          'Indian foot'
    %
    %   Examples
    %   --------
    %   % The following return 'foot'
    %   validateLengthUnit('foot')
    %   validateLengthUnit('feet')
    %   validateLengthUnit('international feet')
    %
    %   % The following return 'kilometer'
    %   validateLengthUnit('kilometer')
    %   validateLengthUnit('km')
    %   validateLengthUnit('kilometre')
    %   validateLengthUnit('kilometers')
    %   validateLengthUnit('kilometres')
    %
    %   % The following return 'U.S. survey foot'
    %   validateLengthUnit('U.S. survey foot')
    %   validateLengthUnit('US survey ft')
    %   validateLengthUnit('sf')
    %   validateLengthUnit('U. S. survey feet')
    %   validateLengthUnit('u s survey foot')
    %
    %   % In the following, a non-char input to validateLengthUnit results in
    %   % an error message referencing a function name, 'FOO', a variable name,
    %   % 'UNIT', and an argument number, 5.
    %   validateLengthUnit(17,'FOO','UNIT',5)
    %
    %   See also UNITSRATIO.
    
    % Copyright 2011-2017 The MathWorks, Inc.
    '''
    # Save a copy of the original UNIT, to be used in case of error.
    originalString = unit
    
    # Preprocess the UNIT:
    #   1. Remove periods and apostrophes.
    #   2. Remove spaces.
    #   3. Convert to lower case.
    unit = unit.replace(".", "")
    unit = unit.replace("'", "")
    unit = unit.replace(" ", "")
    unit = unit.lower()
    
    # Synonyms: Values of UNIT found in the first column corresponded to
    # standard names in the second column.  However, some names in the
    # second column (e.g., 'ussurveyfoot') are not fully in standard form,
    # because they have been contracted and converted to lower case (to match
    # the processing applied to UNIT).  Such names are fully standardized in
    # a subsequent filtering step.
    synonyms = np.array([
        [                 'm',  'meter'],
        [            'meters',  'meter'],
        [             'metre',  'meter'],
        [            'metres',  'meter'],
        [           'microns',  'micron'],
        [                'cm',  'centimeter'],
        [       'centimeters',  'centimeter'],
        [        'centimetre',  'centimeter'],
        [       'centimetres',  'centimeter'],
        [                'mm',  'millimeter'],
        [       'millimeters',  'millimeter'],
        [        'millimetre',  'millimeter'],
        [       'millimetres',  'millimeter'],
        [                'km',  'kilometer'],
        [        'kilometers',  'kilometer'],
        [         'kilometre',  'kilometer'],
        [        'kilometres',  'kilometer'],
        [                'nm',  'nauticalmile'],
        [            'nautmi',  'nauticalmile'],
        [     'nauticalmiles',  'nauticalmile'],
        [                'in',  'inch'],
        [            'inches',  'inch'],
        [                'yd',  'yard'],
        [               'yds',  'yard'],
        [             'yards',  'yard'],
        [                'ft',  'foot'],
        [   'internationalft',  'foot'],
        [ 'internationalfoot',  'foot'],
        [              'feet',  'foot'],
        [ 'internationalfeet',  'foot'],
        [                'mi',  'mile'],
        [             'miles',  'mile'],
        [ 'internationalmile',  'mile'],
        ['internationalmiles',  'mile'],
        [                'sf',  'ussurveyfoot'],
        [        'surveyfoot',  'ussurveyfoot'],
        [          'surveyft',  'ussurveyfoot'],
        [        'ussurveyft',  'ussurveyfoot'],
        [      'ussurveyfoot',  'ussurveyfoot'],
        [        'surveyfeet',  'ussurveyfoot'],
        [      'ussurveyfeet',  'ussurveyfoot'],
        [                'sm',  'ussurveymile'],
        [        'surveymile',  'ussurveymile'],
        [       'surveymiles',  'ussurveymile'],
        [       'statutemile',  'ussurveymile'],
        [      'statutemiles',  'ussurveymile'],
        [     'ussurveymiles',  'ussurveymile'],
        [       'clarkesfoot',  'clarkesfoot'],
        [  'germanlegalmeter',  'germanlegalmetre'],
        [        'indianfoot',  'indianfoot'],
    ])

    if ~np.any(strcmp(unit, synonyms[:, 1])):
        # UNIT doesn't match any of the standard length unit names in the
        # second column of the synonyms array; see if it matches any of the
        # synonyms in the first column. If a match is found, replace the
        # synonym with the standard name.
        index = strcmp(unit, synonyms[:, 0]);
        if any(index):
            # The first column is known to be free of duplicates, so we
            # can assume that exactly one element in index is true.
            # Select the corresponding unit from the second column.
            unit = synonyms[index, 1]
        else:
            raise ValueError(f"validateLengthUnit: unsupported unit: {originalString}")
    
    # Filter contracted names: UNIT now corresponds to a standard name, except
    # that missing spaces, missing periods, and capitalization may need to be
    # restored.
    contractedNames = np.array([
        ['nauticalmile',     'nautical mile'],
        ['ussurveyfoot',     'U.S. survey foot'],
        ['ussurveymile',     'U.S. survey mile'],
        ['clarkesfoot',      'Clarke''s foot'],
        ['germanlegalmetre', 'German legal metre'],
        ['indianfoot',       'Indian foot'],
    ])
    index = strcmp(unit, contractedNames[:, 0])
    if np.any(index):
        # There's a (unique) match with an element of the first column of
        # contractedNames; replace unit with the corresponding element from the
        # second column (the uncontracted form).
        unit = contractedNames[index, 1]
    
    return unit

def unitnamemod(unitNameDisp):
    '''
    %unitnamemod Updates the unit name for the units ratio function 
    %   [UNITNAME,FAC] = unitnamemod(UNITNAME) converts UNITNAME to the
    %   appropriate unit name for use by unitratio, as well as an additional
    %   conversion factor FAC in the case of kft. 
    %   
    %   Unit name changes:
    %      Input Name     |     Output Name     |     Conversion Factor F
    %      kft                  ft                    1e3
    %      nmi                  nm                    1
    %
    %   Allowed units are the same as unitsratio with the added kft option.
    '''
    unitNameToConvert = ['nmi', 'kft']
    unitNameConversion = ['nm', 'ft']
    facConversion = [1, 1e3]
    
    tf = [x.lower() == unitNameDisp.lower() for x in unitNameToConvert]
    
    if any(tf):
        unitNameOut = unitNameConversion[tf]
        fac = facConversion[tf]
    else:
        unitNameOut = unitNameDisp
        fac = 1

    return unitNameOut, unitNameDisp, fac
    
def rhaplot(
    userSetMaxLims,
    vcp,
    vcpangles,
    rngunit,
    rngunitDisp,
    rngfac,
    htunit,
    htunitDisp,
    htfac,
    spow,
    Ns,
    expo,
    anht,
    fc,
    ec,
    hAxes,
    addDataTips,
    rmax,
    hmax,
):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %rhaplot Range-height-angle plot (Blake chart)
    %   rhaplot(userSetMaxLims,vcp,rngunit,rngunitDisp,rngfac,htunit, ...
    %   htunitDisp,htfac,spow,Ns,expo,anht,fc,ec,hAxes,addDataTips)
    %
    %   userSetMaxLims = Logical indicating whether user set rmax and hmax
    %   vcp            = Vertical coverage contour (in rngunit)
    %   vcpangles      = Vertical coverage angles (in degrees)
    %   rngunit        = Range unit of x-axis and vcp used for conversions
    %   rngunitDisp    = Range unit to display in plot
    %   rngfac         = Additional conversion factor to account for kft
    %   htunit         = Height unit of y-axis and anht used for conversions
    %   htunitDisp     = Height unit to display in plot
    %   htfac          = Additional conversion factor to account for kft
    %   spow           = Scale power between 0 and 1
    %   Ns             = Surface refractivity in N-units
    %   expo           = Refraction exponent in 1/km
    %   anht           = Antenna height (in htunit)
    %   fc             = Face (patch) color, N x 3, where N can be 1 or the
    %                    number of columns in vcp
    %   ec             = Edge (contour) color, N x 3, where N can be 1 or the
    %                    number of columns in vcp
    %   hAxes          = Plot axes
    %   addDataTips    = Logical indicating whether data tips should be added.
    %                    Turn off to increase speed.
    %
    %   See the calling function for additional information.
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also blakechart, radarDesigner.
    
    %   Copyright 2021-2024 The MathWorks, Inc.
    '''

    # Nested helper functions
    # --------------------------
    def getheightXY():
        hts = get_heightsamples(hmax, E, anht) * unitsratio('km', htunit)
        srange = slantrange(hts, theta_0, Ns, expo, anhtkm)
        htX, htY = rngth2xy(srange, theta_1, rmaxkm, spow, E)
        
        return htX, htY, hts
    
    def getboundingxy(bdhtx, bdhty, srGeRmax):
        # Solves for the xy coordinate of the boundary of the chart
        if srGeRmax:
            # When range corresponding to hmax exceeds rmax
            # Solve for point where loci for const hmax intersects with const rmax
            if (hmaxkm == rmaxkm): # fzero fails to solve in this case
                intrth0 = pi / 2
            else:
                # make sure answer is between 0 and pi/2
                thsol = brentq(lambda th: slantrange(hmaxkm, th, Ns, expo, anhtkm) - rmaxkm, 0, pi / 2)
                intrth0 = rem(abs(thsol), pi/2)

            intrth1 = atan(E * tan(intrth0)) # Angle on chart
            _, yintr = rngth2xy(rmaxkm, intrth1, rmaxkm, spow, E)
            bdx = [cos(linspace(0, intrth0)), bdhtx[bdhty > yintr],  0]
            bdy = [E * sin(linspace(0, intrth0)), bdhty[bdhty > yintr], 0]
        else:
            # Solve for angle where the const range curve meets y = 1 line
            intrth0 = brentq(lambda th: (E * sin(th) - 1), 0, pi / 2)
            crngth = [theta_0[theta_0 < intrth0], intrth0]
            bdx = [cos(crngth), 0, 0]
            bdy = [E * sin(crngth), 1, 0]

        return bdx, bdy, intrth0
    
    def getangleXY(srGeRmax):
        # Returns X,Y coordinates for angle spokes
        # Set index for angle to be marked
        # theta_0 is 181 samples
        if E <= 5:
            idxa = np.hstack([1, list(range(101, 171, 10)), 181]) - 1 # may need to adjust these for 0 indexing
        elif (E > 5 and E <= 25):
            idxa = np.hstack([list(range(1, 101, 10)), list(range(111, 151, 10)), 181]) - 1 # may need to adjust these for 0 indexing
        elif (E > 25 and E <= 50):
            idxa = np.hstack([list(range(1, 101, 10)), 111, 121, 136, 151, 171, 181]) - 1 # may need to adjust these for 0 indexing
        elif (E > 50 and E<= 75):
            idxa = np.hstack([list(range(1, 101, 10)), 111, 121, 136, 151, 181]) - 1 # may need to adjust these for 0 indexing
        else:
            idxa = np.hstack([list(range(1, 11)), list(range(21, 101, 10)), 111, 181]) - 1 # may need to adjust these for 0 indexing
        
        if srGeRmax:
            theta_0s = theta_0[idxa]
            tha0 = theta_0s[theta_0s <= intrth0] # part with const rng boundary
            st = find(theta_0s > intrth0, 1, 'first')
            thb0 = theta_0[idxa[unlevel(st):]]
            thsz = length(tha0) + length(thb0)
            
            angX = np.vstack([[zeros(thsz)], np.hstack([cos(tha0), htX[-1, idxa[unlevel(st):]]])])
            angY = np.vstack([[zeros(thsz)], np.hstack([E * sin(tha0), htY[-1, idxa[unlevel(st):]]])])
        else:
            theta_0s = theta_0[idxa] # Sampled theta_0
            tha0 = theta_0s[theta_0s <= intrth0] # part with const rng boundary
            thb0 = theta_0s[theta_0s > intrth0] # part with y = 1 boundary
            thb1 = atan(E * tan(thb0)) # Convert to angles on chart for this portion
            thsz = length(tha0) + length(thb1)
            
            angX = np.vstack([[zeros(thsz)], np.hstack([cos(tha0), 1 / tan(thb1)])])
            angY = np.vstack([[zeros(thsz)], np.hstack([E * sin(tha0), ones(1, length(thb1))])])

        return angX, angY
    
    def getrangeXY(srGeRmax):
        # Range related calculations
        rngs = get_rangesamples(rmax, E) * unitsratio('km', rngunit)
        rngX = cos(singcolvec(theta_0)) * ((rngs/rmaxkm) ** spow)
        rngY =  E * sin(singcolvec(theta_0)) * ((rngs / rmaxkm) ** spow)
        
        if srGeRmax:
            # Const range
            for col in range(0, len(rngs)):
                chtbdY = interp1(htX[-1, :], htY[-1, :], rngX[:, col],'pchip') # Boundary # -1's were end
                temp = rngY[:, col]
                temp[temp > chtbdY] = chtbdY[temp > chtbdY]
                rngY[:, col] = temp
        else:
            # Const Range
            rngY[rngY > 1] = 1

        return rngX, rngY, rngs
    
    curveColor = '--mw-graphics-borderColor-axes-tertiary'
    
    # Determine max range and height
    rmax, rmaxkm, rmaxPlotkm, hmax, hmaxkm, hmaxPlotkm = rmaxhmax(userSetMaxLims, vcp, vcpangles, rngunit, rngfac, htunit, htfac, Ns, expo, anht, rmax, hmax)

    # Convert to SI units. All computation done in SI units
    vcpkm = singcolvec(vcp * unitsratio('km', rngunit))
    vcpangles = deg2rad(vcpangles)
    anhtkm = anht * unitsratio('km', htunit)
    
    # Compute ellipticity factor
    E = (rmaxkm / hmaxkm) ** spow # Ref 1, eqn 8
    
    # Define el-angle range (true angle)
    theta_0 = get_elangle() # True angle
    theta_1 = true2chart(theta_0, E) # Angle on chart
    
    # Computation begins
    srGeRmax = slantrange(hmaxkm, 0, Ns, expo, anhtkm) >= rmaxkm
    htX, htY, hts = getheightXY()
    bdx, bdy, intrth0 = getboundingxy(htX[-1, :], htY[-1, :], srGeRmax) # -1's were end
    rngX, rngY, rngs = getrangeXY(srGeRmax)
    angX, angY = getangleXY(srGeRmax)
    
    # Set axes
    xlimMax = (rmaxPlotkm / rmaxkm) ** spow
    ylimMax = (hmaxPlotkm / hmaxkm) ** spow

    fig, hAxes = plt.subplots(figsize=(8, 5))
    plt.axis([0, xlimMax, 0, ylimMax])
    # set(hAxes, 'Color', 'none')
    # hold(hAxes, 'on')
    
    # # Plot background patch
    # patchBackgroundColor = '--mw-graphics-backgroundColor-axes-primary'
    # patchBorderColor = '--mw-graphics-borderColor-axes-primary'
    # hBackground = patch(hAxes, bdx, bdy, 'w', 'HitTest', 'off', 'Tag', 'background')
    # matlab.graphics.internal.themes.specifyThemePropertyMappings(hBackground, 'FaceColor', patchBackgroundColor)
    # matlab.graphics.internal.themes.specifyThemePropertyMappings(hBackground, 'EdgeColor', patchBorderColor)
    # hBackground.Annotation.LegendInformation.IconDisplayStyle = 'off'
    
    # Find constant height x, y points
    if srGeRmax:
        htXGtCosTheta0 = np.array(bsxfun("gt", htX, cos(theta_0))).astype(int)
        idxFind = np.pad(singrowvec(matlabsum(htXGtCosTheta0, 2) >= 1), (0, numel(htY) - numel(singrowvec(matlabsum(htXGtCosTheta0, 2) >= 1))))
        rowVec = np.array(list(range(0, numel(htY))))
        
        for row in rowVec[idxFind]: # Each row -> each const ht curve
            # For each height curve
            idxc = find(htXGtCosTheta0[row, :], 1, 'last')
            
            # Interpolation to make the plot smooth at edge
            htX[row, idxc] = cos(theta_0[idxc])
            htY[row, idxc] =  interp1(htX[row, :], htY[row, :], cos(theta_0[idxc]), 'pchip')
            htY[row, :int(idxc - 1)] = nan
    
    # # Plot custom range, height, angle axes
    plotCustomAxes(hAxes, rngX, rngY, htX, htY, angX, angY, curveColor)
    
    # Create angle labels
    tx = angX[1, :]
    ty = angY[1, :]
    angth, angr = cart2pol(tx, ty)
    angr = angr + 0.01 # Position the label slightly outside the boundary
    angx, angy = pol2cart(angth, angr)
    
    # Add angle labels at edge of plot
    anglbl = rad2deg(chart2true(angth, E))
    formats = ['{:.0f}', '{:.1f}']
    numAng = len(angth)
    numFormat = [1] * numAng
    numFormat = [2 if (anglbl[i] < 1 and anglbl[i] != 0) else 1 for i in range(numAng)]
    angStr = [formats[numFormat[i]-1].format(anglbl[i]) + chr(176) for i in range(numAng)]
    hAngText = text(hAxes, angx, angy, angStr, VerticalAlignment = 'baseline', Clipping = 'on')
    fixOverlappingText(hAngText)
    
    # Add angle labels within the boundary
    if (xlimMax < 1) or (ylimMax < 1): # #ok<BDSCI>
        maxR = sqrt(xlimMax ** 2 + ylimMax ** 2)
        angr = maxR - 0.3 # Position label slightly inside the boundary
        angx, angy = pol2cart(angth, angr)
        hAngInnerText = text(hAxes, angx, angy, angStr, VerticalAlignment = 'baseline', Clipping = 'on', BackgroundColor = 'none', Tag = 'AngleLabels')
        # matlab.graphics.internal.themes.specifyThemePropertyMappings(hAngInnerText,'Color',curveColor)
        fixOverlappingText(hAngInnerText)
    
    # Translate range/angle to Cartesian on chart
    _, c = size(vcpkm)
    vcpX = Nan(numel(vcpangles), c)
    vcpY = Nan(numel(vcpangles), c)
    vcpangc = true2chart(vcpangles, E)
    idxPlot = vcpangc >= 0
    for idx in range(c - 1, -1, -1):
        vcpX[idxPlot, idx], vcpY[idxPlot, idx] = rngth2xy(vcpkm[idxPlot, idx], vcpangc[idxPlot], rmaxkm, spow, E)
    
    # Plot pattern patch
    for idx in range(c - 1, -1, -1):
        idxPlot = ~isnan(vcpX[:, idx]) & ~isnan(vcpY[:, idx])
        # hPatch = patch(hAxes, [[0], [vcpX[idxPlot, idx]], [0]], [[0], [vcpY[idxPlot, idx]], [0]], [1, 1, 1], 'EdgeColor', 'none', 'LineStyle', 'none', 'HitTest', 'off', 'FaceAlpha', 0.6, 'Tag', sprintf('vcpPatch-%d', idx))
        # applyPatchColor(hPatch, fc, idx)
    
    # Plot contour lines
    h_vcp = np.empty(c, dtype=object) 
    for idx in range(c - 1, -1, -1):
        # set endpoints to zero to close fill
        vcpY[:, idx][0] = 0
        vcpY[:, idx][-1] = 0
        vcpX[:, idx][0] = 0
        vcpX[:, idx][-1] = 0
        
        h_vcp[idx] = plt.plot(vcpX[:, idx], vcpY[:, idx])#, 'Tag', sprintf('vcpLine-%d', idx))
        hAxes.fill(vcpX[:, idx], vcpY[:, idx], color='C0', alpha=0.5)

        # applyPlotColor(h_vcp[idx], ec, idx)
    
        # # Plot pattern with no interactions. Interactions will be on custom
        # # data tips.
        # h_vcp(idx).PickableParts = 'none'
    
    #----------- X and Y-Axis Ticks and Labels -----------------------------------
    # Range and height axis labels
    xnorm = (rngs / rmaxkm) ** spow
    ynorm = (hts / hmaxkm) ** spow
    if E <= 5:
        tidx = np.hstack([1, 15, list(range(19, len(xnorm), 2)), len(xnorm)]).astype(int) # may need to adjust these for 0 indexing
        tidy = np.hstack([1, 20, 25, list(range(29, len(ynorm), 2)), len(ynorm)]).astype(int) # may need to adjust these for 0 indexing
    elif (E > 5 and E < 50):
        tidx = np.hstack([1, list(range(10, 18, 2)), list(range(19, len(xnorm), 2)), len(xnorm)]).astype(int) # may need to adjust these for 0 indexing
        tidy = np.hstack([1, list(range(10, 18, 2)), list(range(19, len(ynorm), 2)), len(ynorm)]).astype(int) # may need to adjust these for 0 indexing
    else:
        tidx = np.hstack([list(range(1, len(xnorm), 2)), len(xnorm)]).astype(int) # may need to adjust these for 0 indexing
        tidy = np.hstack([list(range(1, len(ynorm), 2)), len(ynorm)]).astype(int) # may need to adjust these for 0 indexing

    tidx = tidx - 1
    tidy = tidy - 1

    tidx, _ = unique(tidx)
    tidy, _ = unique(tidy)
    
    XTicks = xnorm[tidx]
    XTickLabel = rngs[tidx] * (1 / rngfac) * unitsratio(rngunit, 'km')
    XTicks, XTickLabel = fixOverlappingTicks(hAxes, hAngText, XTicks, XTickLabel, 'x')
    
    YTicks = ynorm[tidy]
    YTickLabel = hts[tidy] * (1 / htfac) * unitsratio(htunit, 'km')
    YTicks, YTickLabel = fixOverlappingTicks(hAxes, hAngText, YTicks, YTickLabel, 'y')
    
    # Set X tick marks
    hAxes.set_xticks(XTicks)
    hAxes.set_xticklabels(XTickLabel)
    hAxes.set_xlabel(f"Range ({rngunitDisp})")
    
    # Set Y tick marks
    hAxes.set_yticks(YTicks)
    hAxes.set_yticklabels(YTickLabel)
    hAxes.set_ylabel(f"Height ({htunitDisp})")
    
    # # Add title
    title_text_object = hAxes.set_title("Blakechart")
    pos = title_text_object.get_position()
    hAxes.set_title("Blakechart", x = pos[0], y = pos[1] + 0.04)
    
    return h_vcp

# Validation/parsing helper functions
# -----------------------------------------+
def rmaxhmax(
    userSetMaxLims,
    vcp,
    vcpangles,
    rngunit,
    rngfac,
    htunit,
    htfac,
    Ns,
    expo,
    anht,
    rmax,
    hmax,
):
    # Set rmax/hmax
    if userSetMaxLims:
        # Adjust for potential kft input
        rmaxPlot = rmax * rngfac
        hmaxPlot = hmax * htfac
        
        # Estimate actual maximum range and height
        rmax, hmax = estrmaxhmax(vcp, vcpangles, anht, rngunit, htunit, Ns, expo, rmax, hmax)
    else:
        # Estimate maximum range and height
        rmax, hmax = estrmaxhmax(vcp, vcpangles, anht, rngunit, htunit, Ns, expo)
        rmaxPlot = rmax
        hmaxPlot = hmax
    
    # Convert to km
    rmaxkm = rmax * unitsratio('km', rngunit)
    hmaxkm = hmax * unitsratio('km', htunit)
    rmaxPlotkm = rmaxPlot * unitsratio('km', rngunit)
    hmaxPlotkm = hmaxPlot * unitsratio('km', htunit)
    anhtkm = anht * unitsratio('km',htunit)
    hmaxkm = hmaxkm + anhtkm
    hmaxPlotkm = hmaxPlotkm + anhtkm
    
    # RMAX >= HMAX
    if rmaxkm < hmaxkm:
        rmaxkm = hmaxkm
        
    return rmax, rmaxkm, rmaxPlotkm, hmax, hmaxkm, hmaxPlotkm
    
# Calculation helper functions
# -----------------------------------------+
def slantrange(height, theta, Ns, d, anht):
    '''
    %SLANTRANGE Coordinate on range-height-angle chart
    %   SRANGE = SLANTRANGE(HEIGHT, THETA) Returns propagated range (SRANGE)
    %   for given HEIGHT and ANGLE on a range-height-angle chart
    %   HEIGHT Height from earth surface in km
    %
    %   Ref:
    %   [1] Radio Ray (Radar) Range-Height-Angle Charts, L. V. Blake, NRL,
    %       1968 Eq 4, p 11.
    %   [2] Blake, L.V. "Ray Height Computation for a Continuous Nonlinear
    %       Atmospheric Refractive-Index Profile." RADIO SCIENCE, Vol. 3 (New
    %       Series), No. 1, Jan. 1968.
    %   [3] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
    %       Atmosphere. Washington, DC: US Gov. Print. Off, 1959.
    '''
    
    # All calculation are done in SI units
    if (numel(height) == 1) and (numel(theta) == 1):
        srange = rangeIntegralCRPL(height, anht, theta, Ns, d)
    else:
        hlog = logspace(-12, 2, int(1e4))
        hint = matlabunion(height, hlog)
        y = rangeIntegrandCRPL(hint, anht, theta, Ns, d)
        srangeCumtrapz = cumtrapz(hint, y)
        idxHeight = [int(i) for i in interp1(hint, list(range(0, numel(hint))), height, 'nearest')]
        srange = srangeCumtrapz[idxHeight, :]
        srange[:, 0] = rangeIntegralCRPL(height, anht, theta[0], Ns, d)

    return srange

def rangeIntegrandCRPL(h, ha, theta0, Ns, k):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %rangeIntegrandCRPL Integrand for propagated range using CRPL 
    %   R = rangeIntegrandCRPL(H,HA,THETA0,NS,K) calculates the propagated
    %   range R in km using the CRPL atmospheric model. 
    %
    %   The inputs are as follows: 
    %      H      = Target height (km)
    %      HA     = Antenna height (km)
    %      THETA0 = Elevation angle at radar (rad)
    %      NS     = Surface refractivity (N-units, 0 height)
    %      K      = CRPL decay constant (1/km)
    %      
    %   See the calling function for additional information. 
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also height2range, height2grndrange, el2height.
    
    %   Copyright 2021 The MathWorks, Inc.
    
    %   References
    %   [1] Blake, L.V. "Ray Height Computation for a Continuous Nonlinear
    %       Atmospheric Refractive-Index Profile." RADIO SCIENCE, Vol. 3 (New
    %       Series), No. 1, Jan. 1968, pp. 85-92.
    %   [2] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
    %       Atmosphere. Washington, DC: US Gov. Print. Off, 1959.
    '''

    h = np.array(singcolvec(h))
    theta0 = np.array(singrowvec(theta0))
    
    r0 = Rearth * 1e-3 + ha # Eqn 8 (km)
    Ns = Ns * 1e-6
    rho0 = Ns * exp(-k * ha) # Eqn 7
    n = 1 + rho0 * exp(-k * h)
    u = (1 + rho0) ** 2 * sin(theta0) ** 2 - 2 * rho0 - rho0 ** 2
    v = 2 * rho0 * exp(-k * h) + rho0 ** 2 * exp(-2 * k * h)
    w = 2 * h / r0 + h ** 2 / r0 ** 2

    R = np.absolute(n ** 2 * (1 + h / r0) / sqrt(u + v + w + v * w)) # Eqn 16
    
    return R

def rangeIntegralCRPL(heightsIn, haIn, theta0In, Ns, k):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %rangeIntegrandCRPL Integrand for propagated range using CRPL 
    %   R = rangeIntegrandCRPL(H,HA,THETA0,NS,K) calculates the propagated
    %   range R in km using the CRPL atmospheric model. 
    %
    %   The inputs are as follows: 
    %      H      = Target height (km)
    %      HA     = Antenna height (km)
    %      THETA0 = Elevation angle at radar (rad)
    %      NS     = Surface refractivity (N-units, 0 height)
    %      K      = CRPL decay constant (1/km)
    %      
    %   See the calling function for additional information. 
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also height2range, height2grndrange, el2height.
    
    %   Copyright 2021-2022 The MathWorks, Inc.
    
    %   References
    %   [1] Blake, L.V. "Ray Height Computation for a Continuous Nonlinear
    %       Atmospheric Refractive-Index Profile." RADIO SCIENCE, Vol. 3 (New
    %       Series), No. 1, Jan. 1968, pp. 85-92.
    %   [2] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
    %       Atmosphere. Washington, DC: US Gov. Print. Off, 1959.
    '''
    
    # Expand vectors as necessary
    maxNum, _ = matlabmax([numel(heightsIn), numel(haIn), numel(theta0In)])
    heights = expandVector(heightsIn, maxNum)
    ha = expandVector(haIn, maxNum)
    theta0 = expandVector(theta0In, maxNum)

    if maxNum == 1:
        heights = np.array([heights])
        ha = np.array([ha])
        theta0 = np.array([theta0])
    
    # Only perform integral for values of theta > 0.1 degrees and heights > 1
    # ft (3.048e-04 km). See pg 89 in reference 1.
    idx = np.array(list(range(0, maxNum)))
    idxNear0 = (abs(theta0) < deg2rad(0.1))
    idxIntegrate = ~(idxNear0 & (abs(heights) < 3.048e-04))
    R = zeros(maxNum)
    if np.any(idxIntegrate):
        # Loosen integration error bounds in the case when angle is near 0
        if np.any(idxNear0):
            thisIdx = idxIntegrate & idxNear0
            R[thisIdx] = np.absolute(arrayfun(lambda x: quadgk(lambda h: rangeIntegrandCRPL(h, ha[x], theta0[x], Ns, k), 0, heights[x], AbsTol = 1e-2, RelTol = 0)[0], idx[thisIdx])) # Eqn 2
        
        # Perform normal calculations for larger angles
        if np.any(~idxNear0):
            thisIdx = idxIntegrate & ~idxNear0
            R[thisIdx] = np.absolute(arrayfun(lambda x: quadgk(lambda h: rangeIntegrandCRPL(h, ha[x], theta0[x], Ns, k), 0, heights[x])[0], idx[thisIdx])) # Eqn 2
    
    # Perform approximation outlined in section 5 of reference 1 for values
    # near theta0 = 0 and h = 0
    if np.any(~idxIntegrate):
        R[idx[~idxIntegrate]] = rangeIntegralApprox(heights[~idxIntegrate], ha[~idxIntegrate], theta0[~idxIntegrate], Ns, k) # Eqn 22

    return R
    
def expandVector(xin, maxNum):
    if isscalar(xin):
        xout = xin * ones((1, maxNum))
        # xout = xin(1).*ones(1,maxNum);
    else:
        xout = xin

    return unlevel(xout)
    
def rangeIntegralApprox(h, ha, theta0, Ns, k):
    r0 = Rearth * 1e-3 + ha
    rho0 = Ns * 1e-6 * exp(-k * ha)
    gammaConst = rho0 * k / (1 + rho0) # Eqn 23
    A_R = 1 + rho0 # Eqn 24
    a_R, _ = matlabmax(k * gammaConst - 7 * gammaConst ** 2 + 8 * gammaConst / r0 - 3 / r0 ** 2 + sin(theta0) ** 2 * (10 * gammaConst ** 2 - 8 * gammaConst / r0 - 2 * k * gammaConst + 3 / r0 ** 2), eps()) 
    b_R = 2 * ((2 * gammaConst - 1 / r0) * sin(theta0) ** 2 - gammaConst + 1 / r0)
    c_R = sin(theta0) ** 2
    R = taylorExpansionApprox(h, A_R, a_R, b_R, c_R)
    
    return R
    
def taylorExpansionApprox(h1, A, a, b, c):
    h1 = abs(h1)
    x = (2 * a * h1 - 2 * sqrt(a * c) + 2 * sqrt(a * (a * h1 ** 2 + b * h1 + c))) / (b + 2 * sqrt(a * c))
    F = abs(A / sqrt(a) * log(1 + x))

    return F

def rngth2xy(rng, th1, rmaxkm, pow_in, E):
    # RNGTH2XY Maps range and angle th1 (angle on chart) to x, y coordinates
    #   Ref: Radio Ray (Radar) Range-Height-Angle Charts, L. V. Blake, NRL,
    #   1968 Eq 4, p 11
    parma = (rng / rmaxkm) ** pow_in # Eqn 11, p 13
    parmb = E * parma # Eqn 12, p 13
    
    # Compute x, y coordinates
    xcor = bsxfun("rdivide", parma, sqrt(1 + (tan(th1) / E) ** 2)) # Eqn 9, p 13
    ycor = bsxfun("times", parmb, sqrt(1 - (xcor / parma) ** 2)) # Eqn 10, p 13
    
    return xcor, ycor
    
def true2chart(theta_tr, E):
    # Converts true angles to angles on chart
    # theta_ch = theta_tr if E = 1
    
    # Eq 13, p 13, NRL Report 6650
    
    theta_ch = atan(E * tan(theta_tr))
    return theta_ch
    
def chart2true(theta_ch, E):
    # Converts true angles to angles on chart
    # theta_ch = theta_tr if E = 1

    # Eq 13, p 13, NRL Report 6650
    
    theta_tr = atan(tan(theta_ch) / E)
    return theta_tr
    
# Sample generation helper functions
# -----------------------------------------+
def estrmaxhmax(vcp, vcpang, anht, rngunit, htunit, Ns, expo, rinput = [], hinput = []):
    # Estimate rmax and hmax limits for blakechart based on values of vcp
    # rmax is estimated using the maximum range covered of vcp
    # hmax is estimated using the mean range coverage of vcp
    # This function needs to be optimized to better estimate hmax.

    if np.ndim(vcp) <= 1:
        vcp = np.array([vcp])

    # Set rmax
    vcpm = vcp * unitsratio('m', rngunit)
    rlim1, _ = matlabmax(vcp[-1, :] * cosd(vcpang)) # Estimate 1 of maximum range 
    rlim1 = nextfactor(rlim1)
    _, idxMaxEst = matlabmax(vcp[-1, :] * sind(vcpang)) # Determine the index of the potential maximum height
    _, idxMinEst = matlabmin(vcp[-1, :] * sind(vcpang)) # Determine the index of the potential minimum height

    rEst = [vcpm[-1, i] for i in [idxMaxEst, idxMinEst]]
    angEst = [vcpang[i] for i in [idxMaxEst, idxMinEst]]
    anhts = repmat(anht * unitsratio('m', htunit) * 1e-3, 1, 2) 
    tgthts = (rEst[:] * sind(angEst[:])).T * 1e-3
    tgthts = tgthts - anhts
    rlim2 = abs(rangeIntegralCRPL(tgthts, anhts, deg2rad(angEst).T, Ns, expo)) * 1e3 * unitsratio(rngunit, 'm')
    rlim2 = nextfactor(matlabmax(rlim2[:])[0])

    if rinput != []:
        rmax, _ = matlabmax(np.vstack([[rlim1], [rlim2], [rinput]]))
    else:
        rmax, _ = matlabmax(np.vstack([[rlim1], [rlim2]]))

    if numel(rmax) > 1:
        rmax[isnan(rmax)] = 100e3 * unitsratio(rngunit, 'm')
    elif isnan(rmax):
        rmax[isnan(rmax)] = 100e3 * unitsratio(rngunit, 'm')
    
    # Set hmax
    hEst = range2height(vcpm[-1, idxMaxEst], anht * unitsratio('m', htunit), vcpang[idxMaxEst], method = 'CRPL', Ns = Ns, rexp = expo, maxNumIter = 1, tol = 1e-4) * unitsratio(htunit, 'm')
    hEst = nextfactor(hEst)
    hmax, _ = matlabmax(np.vstack([singcolvec(hEst), singcolvec(hinput)]))
    hmax = np.array(hmax)
    hmax[isnan(hmax)] = 100e3 * unitsratio(htunit, 'm')
    
    return rmax, hmax
    
def nextfactor(x):
    # Returns the next 10 multiple after x
    # Eg nextfactor(13) returns 20, nextfactor(123) returns 200
    ex = fix(log10(x))
    fac = 10 ** ex
    re = rem(x, fac)
    y = (x - re) / fac
    nf = (y + 1) * fac
    return nf
    
def get_elangle():
    # Returns properly sample elevation angles in deg for RHA chart generation
    # Based on code snippet from NRL report 7098
    # Needs optimization
    idx = 0
    elangle = zeros(181)
    stedIdx = np.array([[0, 100], [10, 91]])
    for p in range(0, 2):
        for k in range(stedIdx[p, 0], stedIdx[p, 1]):
            n = k * 10 ** p
            elev = n / 10
            elangle[idx] = deg2rad(elev)
            idx = idx + 1

    return elangle
    
def get_heightsamples(hmax, E, anht):
    # Returns height sample points in range [0 hmax] for constant height curves
    # Sample points for height scale
    expe = fix(log10(hmax))
    
    if E <= 5:
        exps = expe - 3 # 3 decade scale
    elif (E > 5 and E < 50):
        exps = expe - 2
    else:
        exps = expe - 1
    
    ht0 = float(10) ** exps # Starting point
    ht = zeros(9 * (expe-exps) + 1)
    ht[0] = ht0
    for k in range(exps, expe):
        ht[(k - exps) * 9 + np.arange(1, 10)] = np.linspace(2 * 10 ** k, 10 ** (k + 1), 9)

    htx = ht[-1] + 10 ** (k + 1)
    htsamp = unique(np.hstack([ht, np.arange(htx, hmax + 1, 0.5 * 10 ** (k + 1)), hmax]))[0]
    htsamp = htsamp + anht

    return htsamp
    
def get_rangesamples(rmax, E):
    # Returns range sample points in range [0 rmax] for constant height curves
    # Sample points for height scale
    expe = fix(log10(rmax))
    
    if E <= 5:
        exps = expe - 3 # 3 decade scale
    elif (E > 5 and E < 50):
        exps = expe - 2 # 2 decade
    else:
        exps = expe - 1 # 1 decade
    
    rng0 = float(10) ** exps # Starting point
    rng = zeros(9 * (expe - exps) + 1)
    rng[0] = rng0
    for k in range(exps, expe):
        rng[(k - exps) * 9 + np.arange(1, 10)] = np.linspace(2 * 10 ** k, 10 ** (k + 1), 9)

    rngt = rng[-1]  + 0.5 * 10 ** (k + 1)
    rngsamp = unique(np.hstack([rng, np.arange(rngt, rmax + 1, 10 ** (k + 1)), rmax]))[0]

    return rngsamp

def range2height(
    R,
    anht,
    el,
    method = "Curved", # Method
    Re = (4 / 3) * Rearth, # EffectiveEarthRadius
    Ns = 313, # SurfaceRefractivity
    rexp = 0.143859, # RefractionExponent
    maxNumIter = None, # MaxNumIterations
    tol = 1e-6, # Tolerance
):
    '''
    %range2height Convert propagated range to target height
    %   TGTHT = range2height(R,ANHT,EL) returns the target height in meters
    %   assuming a curved Earth model with a 4/3 effective Earth radius, which
    %   is an approximation used for modeling refraction effects in the
    %   troposphere. R is the propagated range between the target and the
    %   sensor in meters. ANHT is the sensor height in meters. EL is the
    %   initial elevation angle of the ray leaving the sensor, also known as
    %   the local elevation angle, in degrees. Heights are referenced to the
    %   ground.
    %
    %   R, ANHT, and EL can either be scalars or matching M-length vectors.
    %
    %   TGTHT = range2height(...,'Method',METHOD) specifies the model for the
    %   calculation as 'Curved' | 'Flat' | 'CRPL'. The default method is
    %   'Curved'.
    %      - 'Curved'
    %         Assumes a curved Earth model with a 4/3 effective Earth radius,
    %         which is an approximation used for modeling refraction effects in
    %         the troposphere. Alternative values for the effective Earth
    %         radius can be set using the 'EffectiveEarthRadius' name-value
    %         argument. This is the default method.
    %      - 'Flat'
    %         Assumes a flat Earth model. In the case of a flat Earth model,
    %         the effective Earth radius is infinite. 
    %      - 'CRPL'
    %         Assumes a curved Earth model where the atmosphere is defined by
    %         the CRPL exponential reference atmosphere with a refractivity of
    %         313 N-units and a refraction exponent (decay constant) of
    %         0.143859 1/km.
    %
    %   TGTHT = range2height(...,'Method','Curved','EffectiveEarthRadius',RE)
    %   specifies the effective Earth radius as a positive scalar in meters.
    %   This input applies only when 'Method' is set to 'Curved' and is
    %   otherwise ignored. The default calculates the effective Earth radius
    %   using a refractivity gradient of -39e-9, which results in approximately
    %   4/3 of the real Earth radius.
    %
    %   TGTHT = range2height(...,'Method','CRPL','SurfaceRefractivity',NS,
    %   'RefractionExponent',REXP) specifies the surface refractivity NS (in
    %   N-Units) as a non-negative scalar and the refraction exponent factor
    %   (decay constant) REXP (in 1/km) as a non-negative scalar for the CRPL
    %   exponential reference atmosphere model. This input is only applicable
    %   when 'Method' is set to 'CRPL' and is otherwise ignored. The default
    %   value of NS is 313 N-Units, and the default value of REXP is 0.143859
    %   1/km.
    %
    %   The atmospheric refraction model is given by
    %        n(h) = 1 + NS*1e-6*exp(-REXP*h),
    %   where h is the height in kilometers and n is the index of refraction.
    %
    %   TGTHT = range2height(...,'MaxNumIterations',MAXITER) specifies the
    %   maximum number of iterations for the CRPL method as a non-negative
    %   scalar integer. This input acts as a safe guard to prevent endless
    %   iterative calculations. This input applies only when 'Method' is set
    %   to 'CRPL' and is otherwise ignored. Defaults to 10.
    %
    %   MAXITER equal to 0 indicates that a quicker but less accurate
    %   non-iterative CRPL calculation should be performed. The non-iterative
    %   calculation was found to have a maximum height error of 0.056388 m
    %   (0.185 ft) at a target height of 30480 m (100,000 ft) and an elevation
    %   angle equal to 0. The height error for the non-iterative method
    %   decreases with decreasing target height and increasing elevation angle.
    %
    %   TGTHT = range2height(...,'Tolerance',TOL) specifies the numerical
    %   tolerance for the CRPL method at which the iterative process is
    %   terminated. This input applies only when 'Method' is set to 'CRPL'
    %   and 'MaxNumIterations' is greater than 0. Otherwise, this input is
    %   ignored. Defaults to 1e-6.
    %
    %   % Examples:
    %
    %   % Example 1: 
    %   %   Determine the target height in meters assuming a 4/3 effective
    %   %   Earth radius given a range of 300 km, a sensor height of 10 meters,
    %   %   and an elevation angle of 0.5 degrees.
    %   R    = 300e3;
    %   anht = 10;
    %   el   = 0.5;
    %   range2height(R,anht,el)
    %
    %   % Example 2: 
    %   % Calculate the target height assuming a flat Earth, free space
    %   % propagation with a curved Earth, a 4/3 effective Earth radius, and a
    %   % CRPL atmospheric model. Assume a range of 200 km and an antenna 
    %   % height of 100 meters. Plot results over a range of elevation angles 
    %   % from 0 to 5 degrees.
    %   R    = 200e3;
    %   anht = 100;
    %   el   = 0:0.1:5;
    %   
    %   % Calculate target height assuming flat Earth
    %   tgthtFlat = range2height(R,anht,el,'Method','Flat');
    %   
    %   % Calculate target height assuming free space propagation with a curved
    %   % Earth
    %   r0          = physconst('EarthRadius'); 
    %   tgthtFS     = range2height(R,anht,el,'Method','Curved', ...
    %       'EffectiveEarthRadius',r0); 
    % 
    %   % Calculate target height assuming 4/3 effective Earth radius using 
    %   % default values of range2height function. 
    %   tgthtEffRad = range2height(R,anht,el);
    %   
    %   % Calculate target height using CRPL atmospheric model. 
    %   tgthtCRPL   = range2height(R,anht,el,'Method','CRPL');
    %   
    %   % Plot results
    %   figure()
    %   plot(el,[tgthtFlat(:) tgthtFS(:) tgthtEffRad(:)], ...
    %       'LineWidth',1.5)
    %   hold on
    %   plot(el,tgthtCRPL,'--','LineWidth',1.5)
    %   grid on 
    %   xlabel('Elevation Angle (deg)')
    %   ylabel('Target Height (m)')
    %   legend('Flat','Free Space','4/3 Earth','CRPL','Location','Best') 
    %   title('Target Height Estimation')
    %
    %   See also height2range, height2grndrange, el2height, height2el,
    %   refractionexp, radarvcd, blakechart.
    
    %   Copyright 2021-2022 The MathWorks, Inc.  
    
    %   Reference
    %   [1] Blake, L.V. "Ray Height Computation for a Continuous Nonlinear
    %       Atmospheric Refractive-Index Profile." RADIO SCIENCE, Vol. 3 (New
    %       Series), No. 1, Jan. 1968, pp. 85-92.
    %   [2] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
    %       Atmosphere. Washington, DC: US Gov. Print. Off., 1959.
    %   [3] Barton, David K. Radar Equations for Modern Radar. Norwood, MA:
    %       Artech House, 2013.
    '''
    isIterative = False
    if maxNumIter is not None:
        isIterative = True
    
    # Calculate height
    match method.lower():
        case 'curved':
            ht = el2height(el, anht, R, 'Curved', Re)
        case 'flat':
            ht = el2height(el, anht, R, 'Flat')
        case _: # CRPL
            # Convert m to km
            anht = anht * 1e-3 # Antenna height (km)
            R = R * 1e-3 # Target height (km)
            
            # Perform CRPL. Height is relative to antenna. 
            if isIterative:
                ht, _, _ = iterativeCRPL(R, anht, el, Ns, rexp, maxNumIter, tol) # Height (km)
            else:
                ht = nonIterativeCRPL(R, anht, el, Ns, rexp) # Height (km)
            
            # Convert km to m
            ht = np.array(ht) * 1e3 # Height (m)
            
    return ht

def iterativeCRPL(Rrow, anhtRow, elRow, Ns, rexp, maxNumIter, tol):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %iterativeCRPL Target height using iterative CRPL method
    %   HT = iterativeCRPL(R,ANHT,EL,NS,REXP,MAXNUMITER,TOL) calculates the
    %   height HT in km using an iterative method assuming a CRPL atmosphere.
    %
    %   The inputs are as follows:
    %      R          = Range (km). M-length row vector.
    %      ANHT       = Antenna height (km). M-length row vector.
    %      EL         = Elevation angle (deg). M-length row vector.
    %      NS         = Refractivity (N-units). Scalar.
    %      REXP       = Decay constant (1/km). Scalar.
    %      MAXNUMITER = Maximum number of iterations. Scalar.
    %      TOL        = Tolerance. Scalar.
    %
    %   See the calling function for additional information.
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also range2height. 
    
    %   Copyright 2021 The MathWorks, Inc.
    
    %   References
    %   [1] Blake, L.V. "Ray Height Computation for a Continuous Nonlinear
    %       Atmospheric Refractive-Index Profile." RADIO SCIENCE, Vol. 3 (New
    %       Series), No. 1, Jan. 1968, pp. 85-92.
    %   [2] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
    %       Atmosphere. Washington, DC: US Gov. Print. Off, 1959.
    '''

    # Setup
    Reff = (4 / 3) * Rearth
    
    # Calculate hmax
    hmax = el2height(elRow, anhtRow * 1e3, Rrow * 1e3, 'Curved', Rearth) * 1e-3 # Eqn 10
    hmax = hmax - anhtRow # Height relative to antenna
    
    # Calculate hmin
    hmin = el2height(elRow, anhtRow * 1e3, Rrow * 1e3, 'Curved', Reff) * 1e-3 # Eqn 10
    hmin = hmin - anhtRow # Height relative to antenna
    
    # hmax ~= hmin ~= 0
    idxNot0 = (hmax > sqrt(eps())) & (hmin > sqrt(eps()))
    if ~np.any(idxNot0):
        # Return if there are no target heights that are not about 0
        h1n = hmin + anhtRow # Height relative to surface 
        numIter = 0
        R1n = Rrow 
        return h1n, numIter, R1n
    
    # Calculate h11
    idxAdjust = abs(hmax - hmin) <= sqrt(eps())
    if np.any(idxAdjust):
        # Adjust hmax/hmin for case where hmax ~ hmin. This typically occurs at
        # angles close to 90 degrees. 
        hmax[idxAdjust] = hmax[idxAdjust] + 0.1
        hmin[idxAdjust] = hmin[idxAdjust] - 0.1
        
    h11 = np.array(hmin)
    h11[idxNot0] = (hmax[idxNot0] + hmin[idxNot0]) / 2 # Eqn 12
    
    # Calculate R11
    elRow = deg2rad(elRow)
    R11 = np.array(Rrow)
    R11[idxNot0] = rangeIntegralCRPL(h11[idxNot0], anhtRow[idxNot0], elRow[idxNot0], Ns, rexp)
    
    # Calculate h12
    h12 = h11
    h12[idxNot0] = h11[idxNot0] + (hmax[idxNot0] - hmin[idxNot0]) / 4 # Eqn 13
    R12 = np.array(Rrow)
    R12[idxNot0] = rangeIntegralCRPL(h12[idxNot0], anhtRow[idxNot0], elRow[idxNot0], Ns, rexp)
    
    # Setup values for first iteration
    R1n = R12
    R1min1 = R12
    R1min2 = R11
    h1min1 = h12
    h1min2 = h11
    h1n = h1min1
    denom = np.array(R1min2 - R1min1)
    
    # Calculate first tolerance
    T = np.array(abs(R12 - Rrow)) # Eqn 15
    idxComplete = np.array(T <= tol)
    
    # Iterate
    numIter = 0
    while (numIter < maxNumIter) and np.any(~idxComplete):
        # Update iteration count
        numIter = numIter + 1
        
        # Calculate new values
        idxIter = ~idxComplete
        denom[idxIter] = R1min2[idxIter] - R1min1[idxIter]
        if abs(denom) < eps(): # Don't let the denominator become zero
            signVal = ones(size(denom))
            signVal[denom < 1] = -1
            denom = eps() * signVal

        h1n[idxIter] = (Rrow[idxIter] - R1min1[idxIter]) * (h1min2[idxIter] - h1min1[idxIter]) / (denom[idxIter]) + h1min1[idxIter] # Eqn 14
        R1n[idxIter] = rangeIntegralCRPL(h1n[idxIter], anhtRow[idxIter], elRow[idxIter], Ns, rexp) # Eqn 16
        
        # Calculate threshold
        T[idxIter] = abs(R1n[idxIter] - Rrow[idxIter]) # Eqn 15
        
        # Update values for next iteration
        R1min2[idxIter] = R1min1[idxIter]
        R1min1[idxIter] = R1n[idxIter]
        h1min2[idxIter] = h1min1[idxIter]
        h1min1[idxIter] = h1n[idxIter]
        
        # Set output for values that satisfy tolerance
        idxComplete[idxIter] = T[idxIter] <= tol 

    h1n = h1n + anhtRow # Height relative to surface

    return h1n, numIter, R1n

def nonIterativeCRPL(Rrow, anhtRow, elRow, Ns, rexp):
    '''
    %This function is for internal use only. It may be removed in the future.
    
    %nonIterativeCRPL Target height using non-iterative CRPL method
    %   HT = nonIterativeCRPL(R,ANHT,EL,NS,REXP) calculates the height HT in km
    %   using a non-iterative method assuming a CRPL atmosphere.
    %
    %   The inputs are as follows:
    %      R          = Range (km). M-length row vector.
    %      ANHT       = Antenna height (km). M-length row vector.
    %      EL         = Elevation angle (deg). M-length row vector.
    %      NS         = Refractivity (N-units). Scalar. 
    %      REXP       = Decay constant (1/km). Scalar. 
    %
    %   See the calling function for additional information.
    %
    %   This is an internal function, and no validation is performed.
    %   Validation should be performed by the calling function or object.
    %
    %   See also range2height. 
    
    %   Copyright 2021 The MathWorks, Inc.
    
    %   References
    %   [1] Blake, L.V. "Ray Height Computation for a Continuous Nonlinear
    %       Atmospheric Refractive-Index Profile." RADIO SCIENCE, Vol. 3 (New
    %       Series), No. 1, Jan. 1968, pp. 85-92.
    %   [2] Bean, B. R., and G. D. Thayer. CRPL Exponential Reference
    %       Atmosphere. Washington, DC: US Gov. Print. Off, 1959.
    '''    
    # Setup 
    Reff = (4 / 3) * Rearth
    
    # Calculate hmin
    hmin = el2height(elRow, anhtRow * 1e3, Rrow * 1e3, 'Curved', Reff) * 1e-3 # Eqn 10
    hmin = hmin - anhtRow # Height relative to antenna
    
    # Calculate Rhmin
    elRow = deg2rad(elRow)
    Rhmin = rangeIntegralCRPL(hmin, anhtRow, elRow, Ns, rexp) # Eqn 2 using hmin
    
    # Calculate delR
    delR = abs(Rrow - Rhmin) # Eqn 43
    
    # Calculate r1, theta1, and rho1
    r0 = Rearth * 1e-3
    r1 = r0 + hmin # Eqn 39
    rho0 = Ns * 1e-6
    theta1 = acos((1 + rho0) * cos(elRow) / ((1 + rho0 * exp(-rexp * hmin)) * (1 + hmin / r0))) # Eqn 40
    rho1 = rho0 * exp(-rexp * hmin) # Eqn 41
    
    # Calculate b1 and c1
    gammaConst = rho1 * rexp / (1 + rho1) # Eqn 23
    b1 = 2 * ((2 * gammaConst - 1 / r1) * sin(theta1) ** 2 - gammaConst + 1 / r1) # Eqn 26
    c1 = sin(theta1) ** 2 # Eqn 27
    
    # Calculate h1 
    h1 = b1 * (delR / (2 * (1 + rho1))) ** 2 + delR * sqrt(c1) / (1 + rho1) + hmin # Eqn 44
    h1 = h1 + anhtRow # Height relative to surface
    
    return h1

def el2height(
    el,
    anht,
    SR,
    model = 'Curved',
    Re = (4 / 3) * Rearth,
):
    '''
    %el2height Convert target initial elevation angle to height
    %   HT = el2height(EL,ANHT,SR) returns the target height in meters. EL is
    %   the target elevation angle in degrees, ANHT is the sensor height in
    %   meters, and SR is the slant range between the target and the sensor in
    %   meters. The computation assumes a curved Earth model with 4/3 effective
    %   Earth radius, which is an approximation used for modeling refraction
    %   effects in the troposphere. Heights are referenced to the ground.
    %
    %   EL, ANHT, and SR can either be scalars or M-length vectors.
    %
    %   HT = el2height(EL,ANHT,SR,MODEL) specifies the Earth model used to
    %   compute the target height as one of 'Flat' | 'Curved', where the
    %   default is 'Curved'.
    %
    %   HT = el2height(EL,ANHT,SR,MODEL,RE) specifies the effective Earth
    %   radius in meters as a positive scalar, RE, where the default is 4/3 of
    %   the Earth radius. The effective Earth radius is an approximation used
    %   for modeling refraction effects in the troposphere. This input is
    %   ignored when you set Model to 'Flat'.
    %
    %   % Example:
    %   %   Determine the target height in meters given an elevation angle of 
    %   %   0.5 degrees, a sensor height of 10 meters, and a slant range of 300
    %   %   km. 
    %   el   = 0.5;
    %   anht = 10;
    %   SR   = 300e3;
    %   el2height(el,anht,SR)
    %
    %   See also height2el, horizonrange, depressionang, grazingang,
    %   effearthradius.
    
    %   Copyright 2020-2022 The MathWorks, Inc.  
    
    %   Reference
    %   [1] Barton, David K. Radar Equations for Modern Radar. Norwood, MA:
    %       Artech House, 2013.
    '''
    
    if model[0] == 'C':
        # Curved Earth calculation
        val = (Re + anht) ** 2 + SR ** 2 + 2 * (Re + anht) * SR * sind(el) 
        ht = sqrt(val) - Re # Equation 7.9 in Barton (m)
    
        # Assume flat Earth calculation for when el ~ 90
        idx90 = 90 - abs(el) < 1
        htFlat = anht + SR * sind(el) # (m)

        if numel(ht) > 1:
            ht[idx90] = htFlat[idx90]
        elif idx90:
            ht = htFlat
    else:
        # Flat Earth calculation
        ht = anht + SR * sind(el) # (m)
        idxNaN = isnan(ht)
        if numel(anht) > 1:
            ht[idxNaN] = anht[idxNaN] # Edge case when R*sind(theta) = inf*0 = NaN
        else:
            ht[idxNaN] = anht # Edge case when R*sind(theta) = inf*0 = NaN
    
    return ht

def cart2pol(x, y, z = None):
    '''
    %CART2POL Transform Cartesian to polar coordinates.
    %   [TH,R] = CART2POL(X,Y) transforms corresponding elements of data stored
    %   in Cartesian coordinates X,Y to polar coordinates (angle TH and radius
    %   R). The arrays X and Y must have compatible sizes. In the simplest
    %   cases, they can be the same size or one can be a scalar. Two inputs
    %   have compatible sizes if, for every dimension, the dimension sizes of
    %   the inputs are either the same or one of them is 1. TH is returned in
    %   radians.
    %
    %   [TH,R,Z] = CART2POL(X,Y,Z) transforms corresponding elements of data
    %   stored in Cartesian coordinates X,Y,Z to cylindrical coordinates (angle
    %   TH, radius R, and height Z). The arrays X,Y, and Z must have compatible
    %   sizes. TH is returned in radians.
    %
    %   Class support for inputs X,Y,Z:
    %      float: double, single
    %
    %   See also CART2SPH, SPH2CART, POL2CART.
    
    %   Copyright 1984-2021 The MathWorks, Inc. 
    '''
    th = atan2(y, x)
    r = hypot(x, y)

    if z is not None:
        return th, r, z
    else:
        return th, r

def pol2cart(th, r, z = None):
    '''
    %POL2CART Transform polar to Cartesian coordinates.
    %   [X,Y] = POL2CART(TH,R) transforms corresponding elements of data stored
    %   in polar coordinates (angle TH, radius R) to Cartesian coordinates X,Y.
    %   The arrays TH and R must have compatible sizes. In the simplest cases,
    %   they can be the same size or one can be a scalar. Two inputs have
    %   compatible sizes if, for every dimension, the dimension sizes of the
    %   inputs are either the same or one of them is 1. TH must be in radians.
    %
    %   [X,Y,Z] = POL2CART(TH,R,Z) transforms corresponding elements of data
    %   stored in cylindrical coordinates (angle TH, radius R, height Z) to
    %   Cartesian coordinates X,Y,Z. The arrays TH, R, and Z must have
    %   compatible sizes.  TH must be in radians.
    %
    %   Class support for inputs TH,R,Z:
    %      float: double, single
    %
    %   See also CART2SPH, CART2POL, SPH2CART.
    
    %   L. Shure, 4-20-92.
    %   Copyright 1984-2021 The MathWorks, Inc. 
    '''
    x = r * cos(th)
    y = r * sin(th)

    if z is not None:
        return x, y, z
    else:
        return x, y
    
# Plotting helper functions
# -----------------------------------------+
def plotCustomAxes(hAxes, rngX, rngY, htX, htY, angX, angY, curveColor):
    # Plot constant height curves
    hHeightAxes = line(hAxes, htX[:-1, :], htY[:-1, :], HandleVis = 'off', HitTest = 'off', Tag = 'ConstantHeightCurves')
    # matlab.graphics.internal.themes.specifyThemePropertyMappings(hHeightAxes, 'Color', curveColor)
    # numH = numel(hHeightAxes)
    # arrayfun(@(idx) set(hHeightAxes(idx).Annotation.LegendInformation, 'IconDisplayStyle', 'off'), 1:numH)
    
    # Plot constant range curves
    hRangeAxes = line(hAxes,rngX, rngY, HandleVis = 'off', HitTest = 'off', Tag = 'ConstantRangeCurves')
    # matlab.graphics.internal.themes.specifyThemePropertyMappings(hRangeAxes, 'Color', curveColor)
    # numH = numel(hRangeAxes)
    # arrayfun(@(idx) set(hRangeAxes(idx).Annotation.LegendInformation, 'IconDisplayStyle', 'off'), 1:numH)
    
    # Plot angle pokes
    hAngleAxes = line(hAxes, angX, angY, HandleVis = 'off', HitTest = 'off', Tag = 'ConstantAngleLines')
    # matlab.graphics.internal.themes.specifyThemePropertyMappings(hAngleAxes, 'Color', curveColor)
    # numH = numel(hAngleAxes)
    # arrayfun(@(idx) set(hAngleAxes(idx).Annotation.LegendInformation, 'IconDisplayStyle', 'off'), 1:numH)
    
    hAngleAxesEnd = line(hAxes, angX[:, -1], angY[:, -1], HandleVis = 'off', HitTest = 'off', Tag = 'ConstantAngleEndPoints')
    # matlab.graphics.internal.themes.specifyThemePropertyMappings(hAngleAxesEnd, 'Color', '--mw-graphics-borderColor-axes-primary')
    # numH = numel(hAngleAxesEnd)
    # arrayfun(@(idx) set(hAngleAxesEnd(idx).Annotation.LegendInformation, 'IconDisplayStyle', 'off'), 1:numH)
    
def fixOverlappingText(hText):
    # Fix overlapping text for angular axes
    numH = numel(hText)
    idxDelete = zeros(numH)
    iiprev = 0
    if numH > 1:
        for ii in range(1, numH):
            x2 = Extent(hText[ii])[0]
            y2 = Extent(hText[ii])[1]
            w2 = Extent(hText[ii])[2]
            h2 = Extent(hText[ii])[3]
            rectPts = np.array([[x2, y2], [x2, (y2 + h2)], [(x2 + w2), (y2 + h2)], [x2, (y2 + h2)]])
            
            # Check for overlapping x
            minX = Extent(hText[iiprev])[0]
            maxX = Extent(hText[iiprev])[0] + Extent(hText[iiprev])[2]
            idxOverlapX = any((rectPts[:, 0] >= minX) & (rectPts[:, 0] <= maxX))
            
            # Check for overlapping y
            minY = Extent(hText[iiprev])[1]
            maxY = Extent(hText[iiprev])[1] + Extent(hText[iiprev])[3]
            idxOverlapY = any((rectPts[:, 1] >= minY) & (rectPts[:, 1] <= maxY))
            
            delCond = idxOverlapX and idxOverlapY
            idxDelete, iiprev = updateDeleteIdx(delCond, idxDelete, iiprev, ii)
        # delete(hText[idxDelete])
    
def fixOverlappingTicks(hAxes, hText, ticks, tickLabel, type_in):
    # Estimate text dimensions
    match type_in:
        case 'x':
            # Estimate tick text width
            tickLength = Extent(hText[0])[2]
            numChar = numel(hText[0].get_text().replace("'", "").split(","))
            tickLength = tickLength / numChar
            
            # Estimate tick extent in horizontal direction
            formatType = hAxes.ticklabel_format(axis = "x")
            xTickLabelsChar = arrayfun(lambda x: num2str(x, formatType), tickLabel, UniformOutput = False)
            numChars = cellfun(lambda x: numel(x) + 1, xTickLabelsChar) # Add extra character for spacing
            halfTickLength = (tickLength * singcolvec(numChars)) / 2
        case _:
            # Estimate tick text height
            tickLength = Extent(hText[0])[3]
            halfTickLength = (tickLength) / 2
    
    # Calculate text extent
    textExtent = np.hstack([(singcolvec(ticks) - halfTickLength) + 1, (singcolvec(ticks) + halfTickLength) - 1]) # Min/max y2
    
    # Identify overlapping tick marks
    numTicks = numel(ticks)
    idxDelete = ones(numTicks).astype(int) # changed from zeros
    iiprev = 1
    if numTicks > 1:
        for ii in range(1, numTicks):
            # delCond = any((textExtent[ii, :] >= textExtent[iiprev, 0]) & (textExtent[ii, :] <= textExtent[iiprev, 1]))
            delCond = True
            idxDelete, iiprev = updateDeleteIdx(delCond, idxDelete, iiprev, ii)

    idxDelete = [bool(i) for i in idxDelete]
    # Delete overlapping tick marks
    ticks = ticks[idxDelete]
    tickLabel = tickLabel[idxDelete]
    
    return ticks, tickLabel
    
def updateDeleteIdx(delCond, idxDelete, iiprev, ii):
    # Update indices to delete
    if delCond:
        # Add index to list of those to be deleted. Do not update previous
        # index.
        idxDelete[ii] = True
    else:
        # Update previous index
        iiprev = ii

    return idxDelete, iiprev
    
# function applyPatchColor(hPatch,fc,idx)
#     if iscell(fc)
#         % Custom colors 
#         hPatch.FaceColor = fc{idx}; 
#     else
#         % Default colors
#         thisColorIndex = fc(idx); 
#         patchColor = sprintf('--mw-graphics-colorOrder-%d-primary',thisColorIndex); 
#         matlab.graphics.internal.themes.specifyThemePropertyMappings(...
#             hPatch,'FaceColor',patchColor);
    
#     end
#     end
    
# function applyPlotColor(hPlot,ec,idx)
#     if iscell(ec)
#         % Custom colors
#         hPlot.Color = ec{idx}; 
#     else
#         % Default colors
#         thisColorIndex = ec(idx);
#         lineColor = sprintf('--mw-graphics-colorOrder-%d-primary',thisColorIndex); 
#         matlab.graphics.internal.themes.specifyThemePropertyMappings(...
#             hPlot,'Color',lineColor);
#     end
#     end