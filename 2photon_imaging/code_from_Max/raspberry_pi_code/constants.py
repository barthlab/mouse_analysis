# camera_preview.py
CAMERA_RESOLUTION = (1024, 768)

# hardware pins
LICKPORT_PIN = 21
WATER_SOLENOID_PIN = 8
AIRPUFF_SOLENOID_PIN = 10
FAKE_SOLENOID_PIN = 12
ENCODER_A_PIN = 24
ENCODER_B_PIN = 22
AIRPUFF_TTL_PULSE = 11
VIDEO_TTL_PULSE = 36
BUZZER_PIN = 26
OCTOPUS_PIN = (37, 35, 33, 31, 29, 32, 38, 40)  # 12:00, 1:30, 3:00, 4:30, 6:00, 7:30, 9:00, 10:30

# data file prefixes
TTL_PREFIX = "ttl_data_"
PUFF_PREFIX = "puff_data_"
DIST_PREFIX = "distance_data_"
LICK_PREFIX = "lick_data_"
PUPIL_PREFIX = "pupil_data_"
VIDEO_PREFIX = "mouse_video_"
SPEED_PREFIX = "speed_data_"
