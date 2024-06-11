#!/bin/env python3

"""
Encoder test script for raspberry pi
"""

import RPi.GPIO as GPIO
import constants
import time
from operator import xor

total_distance = 0
a_state = True
b_state = True


def A(pin):
    global a_state, b_state, total_distance
    a_state = not a_state
    if xor(a_state, b_state):
        total_distance += 1
    else:
        total_distance -= 1


def B(pin):
    global a_state, b_state, total_distance
    b_state = not b_state
    if xor(a_state, b_state):
        total_distance -= 1
    else:
        total_distance += 1


def setup():
    """Set up all the pins and set their initial values"""
    GPIO.setmode(GPIO.BOARD)
    GPIO.setup(constants.ENCODER_A_PIN, GPIO.IN)
    GPIO.setup(constants.ENCODER_B_PIN, GPIO.IN)

    GPIO.add_event_detect(constants.ENCODER_A_PIN, GPIO.BOTH, callback=A)
    GPIO.add_event_detect(constants.ENCODER_B_PIN, GPIO.BOTH, callback=B)


def main():
    start_time = time.time()
    print(start_time)
    while True:
        time.sleep(0.5)
        print(total_distance/(1250*4), "cycle")


if "__main__" == __name__:
    setup()
    main()
    GPIO.cleanup()
