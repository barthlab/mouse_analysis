#!/bin/env python3

"""
Encoder test script for raspberry pi
"""

import RPi.GPIO as GPIO
import numpy as np
import os
import csv
import os.path as path
import constants
import time


class PositionEncoder:
    def __init__(self, leftPin, rightPin, callback=None):
        self.leftPin = leftPin
        self.rightPin = rightPin
        self.value = 0
        self.state = '00'
        self.direction = None
        self.callback = callback if callback is not None else self.register_history
        self.history = []
        GPIO.setup(self.leftPin, GPIO.IN, pull_up_down=GPIO.PUD_DOWN)
        GPIO.setup(self.rightPin, GPIO.IN, pull_up_down=GPIO.PUD_DOWN)
        GPIO.add_event_detect(self.leftPin, GPIO.BOTH, callback=self.transition_occurred)
        GPIO.add_event_detect(self.rightPin, GPIO.BOTH, callback=self.transition_occurred)

    def register_history(self, value, direction):
        self.history.append([time.time(), value, direction])

    def transition_occurred(self, channel):
        p1 = GPIO.input(self.leftPin)
        p2 = GPIO.input(self.rightPin)
        newState = "{}{}".format(p1, p2)

        if self.state == "00":  # Resting position
            if newState == "01":  # Turned right 1
                self.direction = "R"
            elif newState == "10":  # Turned left 1
                self.direction = "L"

        elif self.state == "01":  # R1 or L3 position
            if newState == "11":  # Turned right 1
                self.direction = "R"
            elif newState == "00":  # Turned left 1
                if self.direction == "L":
                    self.value = self.value - 1
                    if self.callback is not None:
                        self.callback(self.value, self.direction)

        elif self.state == "10":  # R3 or L1
            if newState == "11":  # Turned left 1
                self.direction = "L"
            elif newState == "00":  # Turned right 1
                if self.direction == "R":
                    self.value = self.value + 1
                    if self.callback is not None:
                        self.callback(self.value, self.direction)

        else:  # self.state == "11"
            if newState == "01":  # Turned left 1
                self.direction = "L"
            elif newState == "10":  # Turned right 1
                self.direction = "R"
            elif newState == "00":  # Skipped an intermediate 01 or 10 state, but if we know direction then a turn is complete
                if self.direction == "L":
                    self.value = self.value - 1
                    if self.callback is not None:
                        self.callback(self.value, self.direction)
                elif self.direction == "R":
                    self.value = self.value + 1
                    if self.callback is not None:
                        self.callback(self.value, self.direction)

        self.state = newState

    def getValue(self):
        return self.value


class DisplacementEncoder:
    def __init__(self, leftPin, rightPin, callback=None):
        self.leftPin = leftPin
        self.rightPin = rightPin
        self.value = 0
        self.callback = callback if callback is not None else self.register_history
        self.history = []
        GPIO.setup(self.leftPin, GPIO.IN, pull_up_down=GPIO.PUD_DOWN)
        GPIO.setup(self.rightPin, GPIO.IN, pull_up_down=GPIO.PUD_DOWN)
        GPIO.add_event_detect(self.leftPin, GPIO.BOTH, callback=self.transition_occurred)
        GPIO.add_event_detect(self.rightPin, GPIO.BOTH, callback=self.transition_occurred)

    def register_history(self, value):
        self.history.append([time.time(), value])

    def transition_occurred(self, channel):
        self.value += 1
        if self.callback is not None:
            self.callback(self.value)

    def getValue(self):
        return self.value


def valueChanged(value, direction):
    print("* New value: {}, Direction: {}".format(value, direction))


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
    e1 = PositionEncoder(constants.ENCODER_A_PIN, constants.ENCODER_B_PIN)
    # e1 = Encoder(constants.ENCODER_A_PIN, constants.ENCODER_B_PIN, valueChanged)

    save_dir = path.join(path.dirname(__file__), '..', 'data')
    exp_name = input("Experiment ID: ")
    distance_writer = CSVFile(path.join(save_dir, f"{constants.DIST_PREFIX}{exp_name}.csv"),
                              ['time', 'position', 'direction'])

    for _ in range(10):
        time.sleep(5)
        print("Value is {}".format(e1.getValue()))
    for data in e1.history:
        distance_writer.addrow(data)
    print(e1.history)


if "__main__" == __name__:
    GPIO.setmode(GPIO.BOARD)
    main()
    GPIO.cleanup()
