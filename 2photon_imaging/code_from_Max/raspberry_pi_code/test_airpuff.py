import RPi.GPIO as GPIO
import constants
import time

fake_id = constants.FAKE_SOLENOID_PIN
puff_id = constants.AIRPUFF_SOLENOID_PIN
cycle_ids = [37, 35, 33, 31, 29, 32, 38, 40]


def setup():
    GPIO.setmode(GPIO.BOARD)
    for i in cycle_ids:
        print(i)
        GPIO.setup(i, GPIO.OUT)

    for i in cycle_ids:
        GPIO.output(i, GPIO.HIGH)


def complex_puff(duration, pin_id):
    GPIO.output(pin_id, GPIO.LOW)
    time.sleep(duration)
    GPIO.output(pin_id, GPIO.HIGH)


def main():
    start_time = time.time()
    print(start_time)
    while time.time() - start_time < 60:
        time.sleep(10)

        for i in cycle_ids:
            print(f"Puffing {i}")
            complex_puff(0.5, i)
            time.sleep(2)


if "__main__" == __name__:
    setup()
    main()
    GPIO.cleanup()
