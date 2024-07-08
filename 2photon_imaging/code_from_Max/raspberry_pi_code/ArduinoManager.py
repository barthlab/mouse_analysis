#!/bin/env python3

"""
Air puff stimulator for raspberry pi
"""


import csv
import time
import os
import os.path as path
import numpy as np
import gpiozero
import RPi.GPIO as GPIO
import picamera
import constants
import Modules
import argparse
from position_recorder import PositionEncoder, DisplacementEncoder
from operator import xor


args = argparse.ArgumentParser()
args.add_argument("-M", "-Module", type=str, help=f'The module from Modules.{dir(Modules)}')
args.add_argument("-cam", "--cam", action='store_true', help="enable camera recording (require more disk space)")
args.add_argument("-pos", action='store_true', help='record precise distance on the treadmill')
cfg = args.parse_args()
print(cfg)


save_dir = path.join(path.dirname(__file__), '..', 'data')
print(f"Data saved at: {save_dir}")
octopus_stimulus = True
print(f"Using Octopus Stimuli" if octopus_stimulus else "Using Puff-Blank System")
recording = cfg.cam
if recording:
    print("-"*8, "camera recording... pay attention to disk space", "-"*8)
positioning = cfg.pos
if positioning:
    print("-"*8, "position recording... pay attention to voltage consumption", "-"*8)


lick_times = []  # Time that the mouse licked


def get_time_ms():
    return time.monotonic_ns()/1e6


class PiCameraRecordingContextManager:
    def __init__(self, filename):
        self._filename = filename

    def __enter__(self):
        self._camera = None
        if recording:
            self._camera = picamera.PiCamera()
            self._camera.resolution = constants.CAMERA_RESOLUTION
            self._camera.start_recording(self._filename)
        return self._camera

    def __exit__(self, exc_type, exc_value, exc_tb):
        if recording:
            self._camera.stop_recording()
            self._camera.close()
        return None


def tone_generate(buzzer, hz=400, dur=0.5):
    t, time_slot = 0, 0.5*dur/hz
    while t < dur:
        buzzer.on()
        time.sleep(time_slot)
        buzzer.off()
        time.sleep(time_slot)
        t += 1/hz


def triangle_tone(buzzer, start_frequency, end_frequency):
    for buzzer_f in range(start_frequency, end_frequency):
        slot = 0.5 / buzzer_f
        buzzer.on()
        time.sleep(slot)
        buzzer.off()
        time.sleep(slot)


def start_lick(pin):
    global lick_times
    lick_times.append(["start", get_time_ms()])


def stop_lick(pin):
    global lick_times
    lick_times.append(["stop", get_time_ms()])


def setup():
    """Set up all the pins and set their initial values"""
    GPIO.setmode(GPIO.BOARD)
    GPIO.setup(constants.LICKPORT_PIN, GPIO.IN)
    GPIO.setup(constants.WATER_SOLENOID_PIN, GPIO.OUT)
    if octopus_stimulus:
        for pin_id in constants.OCTOPUS_PIN:
            GPIO.setup(pin_id, GPIO.OUT)
    else:
        GPIO.setup(constants.AIRPUFF_SOLENOID_PIN, GPIO.OUT)
        GPIO.setup(constants.FAKE_SOLENOID_PIN, GPIO.OUT)
    GPIO.setup(constants.AIRPUFF_TTL_PULSE, GPIO.OUT)
    GPIO.setup(constants.VIDEO_TTL_PULSE, GPIO.OUT)

    GPIO.output(constants.WATER_SOLENOID_PIN, GPIO.HIGH)
    if octopus_stimulus:
        for pin_id in constants.OCTOPUS_PIN:
            GPIO.output(pin_id, GPIO.HIGH)
    else:
        GPIO.output(constants.AIRPUFF_SOLENOID_PIN, GPIO.HIGH)
        GPIO.output(constants.FAKE_SOLENOID_PIN, GPIO.HIGH)
    GPIO.output(constants.AIRPUFF_TTL_PULSE, GPIO.LOW)
    GPIO.output(constants.VIDEO_TTL_PULSE, GPIO.LOW)

    button = gpiozero.Button(constants.LICKPORT_PIN)
    button.when_pressed = start_lick
    button.when_released = stop_lick


class CSVFile:
    def __init__(self, file_dir, headers):
        self.file_dir = file_dir
        with open(file_dir, 'w') as f:
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(headers)

    def addrow(self, data):
        with open(self.file_dir, 'a') as f:
            writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(data)


def main():
    global distance
    global lick_times

    module = eval(f"Modules.{cfg.M}()")
    exp_name = input("Experiment ID: ")
    buzzer = gpiozero.Buzzer(constants.BUZZER_PIN)
    ttl_writer = CSVFile(path.join(save_dir, f"{constants.TTL_PREFIX}{exp_name}.csv"),
                         ["type", "time", "direction"] if octopus_stimulus else ["type", "time"])
    dist_writer = CSVFile(path.join(save_dir, f"{constants.DIST_PREFIX}{exp_name}.csv"), ['time', 'position', 'direction'] if positioning else ['time', 'displacement'])
    lick_writer = CSVFile(path.join(save_dir, f"{constants.LICK_PREFIX}{exp_name}.csv"), ['time'])

    if positioning:
        e1 = PositionEncoder(constants.ENCODER_B_PIN, constants.ENCODER_A_PIN)
    else:
        e1 = DisplacementEncoder(constants.ENCODER_B_PIN, constants.ENCODER_A_PIN)

    ttl, sti_on, sti_off, water_on, water_off, sti_direction = -1, -1, -1, -1, -1, -1
    solenoid_flag, water_flag, beep_flag = None, None, None
    with PiCameraRecordingContextManager(path.join(save_dir, f"{constants.VIDEO_PREFIX}{exp_name}.h264")) as camera:
        for _, command in enumerate(module.run()):
            if command == 'ShortPulse':
                ttl = get_time_ms()
                GPIO.output(constants.VIDEO_TTL_PULSE, GPIO.HIGH)
                GPIO.output(constants.VIDEO_TTL_PULSE, GPIO.LOW)
                ttl_writer.addrow(["video", ttl])
                ttl = -1
            elif command == 'CheckCamera':
                if recording:
                    camera.wait_recording(0)
            elif command == 'PuffOn':
                assert not octopus_stimulus
                sti_on = get_time_ms()
                GPIO.output(constants.AIRPUFF_SOLENOID_PIN, GPIO.LOW)
                GPIO.output(constants.AIRPUFF_TTL_PULSE, GPIO.HIGH)
                solenoid_flag = 'Puff'
            elif command == "PuffOff":
                GPIO.output(constants.AIRPUFF_SOLENOID_PIN, GPIO.HIGH)
                GPIO.output(constants.AIRPUFF_TTL_PULSE, GPIO.LOW)
                sti_off = get_time_ms()
                solenoid_flag = 'Puff'
            elif command == 'BlankOn':
                assert not octopus_stimulus
                sti_on = get_time_ms()
                GPIO.output(constants.FAKE_SOLENOID_PIN, GPIO.LOW)
                GPIO.output(constants.AIRPUFF_TTL_PULSE, GPIO.HIGH)
                solenoid_flag = 'Blank'
            elif command == "BlankOff":
                GPIO.output(constants.FAKE_SOLENOID_PIN, GPIO.HIGH)
                GPIO.output(constants.AIRPUFF_TTL_PULSE, GPIO.LOW)
                sti_off = get_time_ms()
                solenoid_flag = 'Blank'
            elif command[0] == "PuffOn":
                assert octopus_stimulus
                sti_on = get_time_ms()
                sti_direction = command[1]
                GPIO.output(constants.OCTOPUS_PIN[command[1]], GPIO.LOW)
                GPIO.output(constants.AIRPUFF_TTL_PULSE, GPIO.HIGH)
                solenoid_flag = 'Puff'
            elif command[0] == "PuffOff":
                GPIO.output(constants.OCTOPUS_PIN[command[1]], GPIO.HIGH)
                GPIO.output(constants.AIRPUFF_TTL_PULSE, GPIO.LOW)
                sti_off = get_time_ms()
                sti_direction = command[1]
                solenoid_flag = 'Puff'
            elif command == 'WaterOn':
                water_on = get_time_ms()
                GPIO.output(constants.WATER_SOLENOID_PIN, GPIO.LOW)
                water_flag = "Water"
            elif command == 'WaterOff':
                GPIO.output(constants.WATER_SOLENOID_PIN, GPIO.HIGH)
                water_off = get_time_ms()
                water_flag = "Water"
            elif command == 'NoWaterOn':
                water_on = get_time_ms()
                GPIO.output(constants.FAKE_SOLENOID_PIN, GPIO.LOW)
                water_flag = "NoWater"
            elif command == 'NoWaterOff':
                GPIO.output(constants.FAKE_SOLENOID_PIN, GPIO.HIGH)
                water_off = get_time_ms()
                water_flag = "NoWater"
            elif command[0] == 'Beep':
                hz, dur = command[1], command[2]
                beep_on = get_time_ms()
                GPIO.output(constants.VIDEO_TTL_PULSE, GPIO.HIGH)
                tone_generate(buzzer, hz, dur)
                GPIO.output(constants.VIDEO_TTL_PULSE, GPIO.LOW)
                beep_off = get_time_ms()
                beep_flag = "Beep"
            elif command[0] == "Tone":
                start_hz, end_hz = command[1], command[2]
                beep_on = get_time_ms()
                GPIO.output(constants.VIDEO_TTL_PULSE, GPIO.HIGH)
                triangle_tone(buzzer, start_hz, end_hz)
                GPIO.output(constants.VIDEO_TTL_PULSE, GPIO.LOW)
                beep_off = get_time_ms()
                beep_flag = "Tone"
            elif command[0] == "NoBeep":
                dur = command[1]
                beep_on = get_time_ms()
                GPIO.output(constants.VIDEO_TTL_PULSE, GPIO.HIGH)
                time.sleep(dur)
                GPIO.output(constants.VIDEO_TTL_PULSE, GPIO.LOW)
                beep_off = get_time_ms()
                beep_flag = "NoBeep"
            elif command == 'RegisterSti':
                if beep_flag is not None:
                    if beep_on > -1:
                        ttl_writer.addrow([beep_flag+"On", beep_on])
                    if beep_off > -1:
                        ttl_writer.addrow([beep_flag+"Off", beep_off])
                if solenoid_flag is not None:
                    if sti_on > -1:
                        ttl_writer.addrow([solenoid_flag+"On", sti_on, sti_direction] if octopus_stimulus
                                          else [solenoid_flag+"On", sti_on])
                    if sti_off > -1:
                        ttl_writer.addrow([solenoid_flag+"Off", sti_off, sti_direction] if octopus_stimulus
                                          else [solenoid_flag+"Off", sti_off])
                if water_flag is not None:
                    if water_on > -1:
                        ttl_writer.addrow([water_flag+"On", water_on])
                    if water_off > -1:
                        ttl_writer.addrow([water_flag+"Off", water_off])
                ttl, sti_on, sti_off, water_on, water_off, sti_direction = -1, -1, -1, -1, -1, -1
                solenoid_flag, water_flag, beep_flag = None, None, None
            elif command == 'RegisterBehavior':

                data_sheet = e1.history
                for data in data_sheet:
                    dist_writer.addrow(data)
                e1.history = e1.history[len(data_sheet):]

                prev_lick_times, lick_times = lick_times, []
                for data in prev_lick_times:
                    lick_writer.addrow(data)
            else:
                print(command)
                raise NotImplementedError


if "__main__" == __name__:
    setup()
    main()
    GPIO.cleanup()
