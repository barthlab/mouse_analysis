from gpiozero import Buzzer
from time import sleep


def triangle_sound(start_frequency, end_frequency, duration):
	num_buzzer = int(2*duration/(1/start_frequency + 1/end_frequency))
	for buzzer_id in range(num_buzzer+1):
		ratio = (num_buzzer - buzzer_id)/num_buzzer
		slot = 0.5 * (ratio/start_frequency + (1-ratio)/end_frequency)
		buzzer.on()
		sleep(slot)
		buzzer.off()
		sleep(slot)


buzzer = Buzzer(26)
triangle_sound(200, 1000, 3)
sleep(3)
triangle_sound(200, 500, 1.5)
sleep(3)
triangle_sound(200, 1000, 1.5)
sleep(3)
triangle_sound(200, 500, 3)

