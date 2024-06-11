import time
import random
import numpy as np
from colorist import Color


total_distance = 0
sti_dur = 0.5
print(f"Stimulus duration: {sti_dur} seconds")
water_dur = 0.01
print(f"Water drip duration: {water_dur} seconds")
beep_offset = 0.1
print(f"Duration between beep and puff: {beep_offset} seconds")


def vis(flag, c=None):
    if c == "R":
        print(f"{Color.RED}{flag}{Color.OFF}")
    elif c == "B":
        print(f"{Color.BLUE}{flag}{Color.OFF}")
    elif c == "G":
        print(f"{Color.GREEN}{flag}{Color.OFF}")
    else:
        print(f"{flag}")


class TF:
    def __init__(self):
        # Delay between running script and first trial start
        self.initial_delay = 100  # seconds
        # Delay between end of final trial and program termination
        self.final_delay = 100  # seconds
        # Air puff duration
        self.air_time = 0.5  # seconds
        # Water drip duration
        self.water_time = 0.01  # seconds
        # Time between end of air puff and beginning of water release
        self.air_puff_to_water_release_time = 0.8  # seconds
        # Time between last puff of one train and first puff of other train
        self.inter_puff_delay = 19.5  # seconds
        # Number of puffs per train
        self.num_puffs_in_train = 20
        # Number of trains in a trial
        self.num_trains_in_trial = 1
        # Probability of receiving water
        self.water_prob = 50  # %

    def run(self):
        yield 'ShortPulse'
        time.sleep(self.initial_delay)
        for train_counter in range(self.num_trains_in_trial):
            for puff_counter in range(self.num_puffs_in_train):
                yield 'CheckCamera'
                if random.random() * 100 < self.water_prob:
                    vis("True", 'B')
                    yield "PuffOn"
                    time.sleep(self.air_time)
                    yield "PuffOff"
                else:
                    vis("Fake", 'R')
                    yield "BlankOn"
                    time.sleep(self.air_time)
                    yield "BlankOff"

                # time.sleep(self.air_puff_to_water_release_time)
                #
                # yield "WaterOn"
                # time.sleep(self.water_time)
                # yield "WaterOff"

                time.sleep(self.inter_puff_delay)
                yield 'RegisterSti'
                yield 'RegisterBehavior'
            time.sleep(self.final_delay)
            yield 'RegisterBehavior'

    
class ASSA:
    def __init__(self, pos_prob=0.5, neg_prob=0.0, A2S=True):
        self.initial_delay = 20
        print(f"Delay before first trial starts: {self.initial_delay} seconds")
        self.final_delay = 20
        print(f"Delay after last trial ends: {self.final_delay} seconds")
        self.inter_trial_delay_range = (35, 55)
        print(f"Duration between two trial varied between: {self.inter_trial_delay_range} [seconds]")
        self.time2halt = 10*60 - self.inter_trial_delay_range[1]
        print(f"Time to load final delay: {self.time2halt} seconds")
        self.beep_offset = 1
        print(f"Duration between Beep and Puff: {self.beep_offset} seconds")
        self.pos_prob, self.neg_prob = pos_prob, neg_prob
        print(f"Positive Prediction Probability: {self.pos_prob}, Negative Prediction Probability: {self.neg_prob}")
        self.a2s = A2S
        if self.a2s:
            self.cause_true, self.cause_false = "Beep", "NoBeep"
            self.effect_true, self.effect_false = "Puff", "Blank"
            print(f"Beep is given first, and Puff is given later")
        else:
            self.cause_true, self.cause_false = "Puff", "Blank"
            self.effect_true, self.effect_false = "Beep", "NoBeep"
            print(f"Puff is given first, and Beep is given later")

    def run(self):
        yield 'ShortPulse'
        time.sleep(self.initial_delay)

        t_list = []
        t = self.initial_delay
        while t+self.final_delay < self.time2halt:
            random_delay = np.random.randint(*self.inter_trial_delay_range)
            t += random_delay + self.beep_offset
            t_list.append(random_delay)
        n_trial = len(t_list)
        sti_list = np.zeros(n_trial)
        error = np.random.choice(n_trial, int(n_trial*(self.pos_prob+self.neg_prob)), replace=False)
        pos_error = np.random.choice(error, int(n_trial*self.pos_prob), replace=False)
        sti_list[pos_error] += 1
        sti_list[error] += 1

        for trial_id in range(n_trial):
            yield 'CheckCamera'
            self.vis(f"Random Delay {t_list[trial_id]}s")
            time.sleep(t_list[trial_id])
            if sti_list[trial_id] == 0:
                #  Beep + Puff
                self.vis(self.cause_true)
                yield self.cause_true
                time.sleep(self.beep_offset)
                self.vis(self.effect_true)
                yield self.effect_true
            elif sti_list[trial_id] == 1:
                #  Negative Prediction Trial
                self.vis(self.cause_true)
                yield self.cause_true
                time.sleep(self.beep_offset)
                self.vis(self.effect_false)
                yield self.effect_false
            elif sti_list[trial_id] == 2:
                #  Positive Prediction Trial
                self.vis(self.cause_false)
                yield self.cause_false
                time.sleep(self.beep_offset)
                self.vis(self.effect_true)
                yield self.effect_true
            else:
                raise NotImplementedError

            yield 'RegisterSti'
            yield 'RegisterBehavior'
        time.sleep(self.final_delay)
        yield 'RegisterBehavior'

    def vis(self, flag):
        if flag == 'Puff':
            print(f"{Color.RED}{flag}{Color.OFF}")
        elif flag == 'Blank':
            print(f"{Color.GREEN}{flag}{Color.OFF}")
        elif flag == 'Beep':
            print(f"{Color.YELLOW}{flag}{Color.OFF}")
        elif flag == 'NoBeep':
            print(f"{Color.BLUE}{flag}{Color.OFF}")
        else:
            print(f"{flag}")


#  Cued Calibration
class CCBP(ASSA):
    def __init__(self):
        print("Module [Cued Calibration Beep-Puff] loading...")
        super().__init__(pos_prob=0., neg_prob=0., A2S=True)


class CCPB(ASSA):
    def __init__(self):
        print("Module [Cued Calibration Puff-Beep] loading...")
        super().__init__(pos_prob=0., neg_prob=0., A2S=False)


#  UnCued Calibration
class UCCBP(ASSA):
    def __init__(self):
        print("Module [UnCued Calibration Beep-Puff] loading...")
        super().__init__(pos_prob=1., neg_prob=0., A2S=True)


class UCCPB(ASSA):
    def __init__(self):
        print("Module [UnCued Calibration Puff-Beep] loading...")
        super().__init__(pos_prob=1., neg_prob=0., A2S=False)


#  Random Test
class RTBP(ASSA):
    def __init__(self):
        print("Module [Random Test Beep-Puff] loading...")
        super().__init__(pos_prob=0.5, neg_prob=0., A2S=True)


class RTPB(ASSA):
    def __init__(self):
        print("Module [Random Test Puff-Beep] loading...")
        super().__init__(pos_prob=0.5, neg_prob=0., A2S=False)


# Running Mice
class RM:
    def __init__(self):
        self.initial_delay = 20
        print(f"Delay before first trial starts: {self.initial_delay} seconds")
        self.final_delay = 20
        print(f"Delay after last trial ends: {self.final_delay} seconds")
        self.inter_time = 10*60 - self.initial_delay - self.final_delay
        print(f"Experiment duration: {self.inter_time} seconds")
        self.cycles = 20
        print(f"Mice need to run for {self.cycles} cycles")
        self.safe_area = 2
        print(f"Wheel is uniformly dived into {self.safe_area*2} area, odd areas are safe, start from 0")
        self.short_break = 0.05
        print(f"Detection Interval is {self.short_break*1000} ms")
        self.start_time = None
        self.rolling_speed = self.cycles/self.inter_time

    def check_safe(self, absolute_distance):
        relative_distance = (absolute_distance/(1250*4)) % 1
        normalize_distance = (relative_distance - (time.time()-self.start_time)*self.rolling_speed) % 1
        return int(np.floor(normalize_distance*(2*self.safe_area))) % 2 == 0

    def run(self):
        global total_distance

        yield 'ShortPulse'
        time.sleep(self.initial_delay)
        yield 'CheckCamera'

        # start signal
        yield "Beep", 400, 0.5
        time.sleep(0.5)
        yield "Beep", 400, 0.5

        self.start_time = time.time()
        puff_flag = False

        while time.time()-self.start_time < self.inter_time:
            time.sleep(self.short_break)
            safe_flag = self.check_safe(total_distance)
            self.vis(total_distance, safe_flag)
            if safe_flag:
                if puff_flag is True:
                    yield "PuffOff"
                    puff_flag = False
                    yield 'RegisterSti'
                    yield 'RegisterBehavior'
            else:
                if puff_flag is False:
                    yield "PuffOn"
                    puff_flag = True
                    yield 'RegisterSti'
                    yield 'RegisterBehavior'
        if puff_flag is True:
            yield "PuffOff"
            puff_flag = False
            yield 'RegisterSti'
            yield 'RegisterBehavior'

        # end signal
        yield "Beep", 400, 0.5
        time.sleep(0.5)
        yield "Beep", 400, 0.5
        time.sleep(0.5)
        yield "Beep", 400, 0.5

        time.sleep(self.final_delay)
        yield 'RegisterBehavior'

    def vis(self, dis, flag=None):
        normalize_distance = dis/(1250*4)
        if flag is True:
            print(f"Distance: {normalize_distance:.2} Zone:{Color.GREEN}{flag}{Color.OFF}")
        elif flag is False:
            print(f"Distance: {normalize_distance:.2} Zone:{Color.RED}{flag}{Color.OFF}")
        else:
            print(f"Distance: {normalize_distance:.2}")


# Induced Prediction
class IP:
    def __init__(self, pos_prob=0.2, neg_prob=0.4):
        self.initial_delay = 60
        print(f"Delay before first trial starts: {self.initial_delay} seconds")
        self.final_delay = 60
        print(f"Delay after last trial ends: {self.final_delay} seconds")
        self.pos_prob, self.neg_prob = pos_prob, neg_prob
        print(f"Positive Prediction Probability: {self.pos_prob}, Negative Prediction Probability: {self.neg_prob}")
        self.inter_trial_delay_range = (15, 25)
        print(f"Duration between two trial varied between: {self.inter_trial_delay_range} [seconds]")
        self.time2halt = 10*60 - self.initial_delay - self.final_delay
        print(f"Testing duration: {self.time2halt} seconds")

    def run(self):
        yield 'ShortPulse'
        time.sleep(self.initial_delay)

        t_list = []
        t = 0
        while t < self.time2halt:
            random_delay = np.random.randint(*self.inter_trial_delay_range)
            t += random_delay + 5
            t_list.append(random_delay)
        n_trial = len(t_list)
        sti_list = np.zeros(n_trial)
        error = np.random.choice(n_trial, int(n_trial * (self.pos_prob + self.neg_prob)), replace=False)
        pos_error = np.random.choice(error, int(n_trial * self.pos_prob), replace=False)
        sti_list[pos_error] += 1
        sti_list[error] += 1

        for trial_id in range(n_trial):
            yield 'CheckCamera'
            vis(f"Random Delay {t_list[trial_id]}s")
            time.sleep(t_list[trial_id])
            if sti_list[trial_id] == 0:
                #  Beep + Puff
                vis("Tone + Puff", "G")
                yield "Tone", 200, 1600
                time.sleep(beep_offset)
                yield "PuffOn"
                time.sleep(sti_dur)
                yield "PuffOff"
            elif sti_list[trial_id] == 1:
                #  Negative Prediction Trial
                vis("Tone + NoPuff", "B")
                random_end = np.random.choice([400, 600])
                yield "Tone", 200, random_end
                time.sleep(beep_offset)
                yield "BlankOn"
                time.sleep(sti_dur)
                yield "BlankOff"
            elif sti_list[trial_id] == 2:
                #  Positive Prediction Trial
                vis("NoTone + Puff", "R")
                yield "NoBeep", beep_offset
                yield "PuffOn"
                time.sleep(sti_dur)
                yield "PuffOff"
            else:
                raise NotImplementedError

            yield 'RegisterSti'
            yield 'RegisterBehavior'
        time.sleep(self.final_delay)
        yield 'RegisterBehavior'


# two direction
class TD:
    def __init__(self):
        self.initial_delay = 60
        print(f"Delay before first trial starts: {self.initial_delay} seconds")
        self.final_delay = 60
        print(f"Delay after last trial ends: {self.final_delay} seconds")
        self.random = False
        if self.random:
            print(f"Stimulus is provide in a pseudo random way")
        else:
            print(f"Stimulus is provide in a cyclic way")
        self.inter_trial_delay_range = (18, 22)
        print(f"Duration between two trial varied between: {self.inter_trial_delay_range} [seconds]")
        self.time2halt = 10 * 60 - self.initial_delay - self.final_delay
        print(f"Testing duration: {self.time2halt} seconds")

    def run(self):
        yield 'ShortPulse'
        time.sleep(self.initial_delay)

        t_list = []
        t = 0
        while t < self.time2halt:
            random_delay = np.random.randint(*self.inter_trial_delay_range)
            t += random_delay + sti_dur
            t_list.append(random_delay)

        n_trial = len(t_list)
        sti_list = np.zeros(n_trial)
        if self.random:
            second_sti = np.random.choice(n_trial, int(n_trial/2), replace=False)
            sti_list[second_sti] += 1
        else:
            second_sti = np.arange(int(np.ceil(n_trial/2)))*2
            sti_list[second_sti] += 1

        for trial_id in range(n_trial):
            yield 'CheckCamera'
            vis(f"Random Delay {t_list[trial_id]}s")
            time.sleep(t_list[trial_id])
            if sti_list[trial_id] == 0:
                #  first sti: Puff
                vis("Vertical", "G")
                yield "PuffOn"
                time.sleep(sti_dur)
                yield "PuffOff"
            elif sti_list[trial_id] == 1:
                #  second sti: Blank
                vis("Horizontal", "B")
                yield "BlankOn"
                time.sleep(sti_dur)
                yield "BlankOff"
            else:
                raise NotImplementedError

            yield 'RegisterSti'
            yield 'RegisterBehavior'
        time.sleep(self.final_delay)
        yield 'RegisterBehavior'


# Eight direction
class ED:
    def __init__(self):
        print("Module [Eight Directions] loading...")
        self.initial_delay = 60
        print(f"Delay before first trial starts: {self.initial_delay} seconds")
        self.final_delay = 30
        print(f"Delay after last trial ends: {self.final_delay} seconds")
        self.random = True
        if self.random:
            print(f"Stimulus is provide in a pseudo random way")
        else:
            print(f"Stimulus is provide in a cyclic way")
        self.inter_trial_delay_range = (8, 12)
        print(f"Duration between two trial varied between: {self.inter_trial_delay_range} [seconds]")

    def run(self):
        yield 'ShortPulse'
        time.sleep(self.initial_delay)

        if self.random:
            directions = list(np.concatenate([np.random.permutation(8) for _ in range(6)]))
        else:
            directions = [i for _ in range(6) for i in range(8)]
        print(directions)

        t_list = []
        t = 0
        for _ in range(6*8):
            random_delay = np.random.randint(*self.inter_trial_delay_range)
            t += random_delay + sti_dur
            t_list.append(random_delay)

        n_trial = len(t_list)
        for trial_id in range(n_trial):
            yield 'CheckCamera'
            vis(f"Random Delay {t_list[trial_id]}s, Direction: {directions[trial_id]}")

            time.sleep(t_list[trial_id])
            vis("Puff", "G" if trial_id % 2 == 0 else "B")
            yield "PuffOn", directions[trial_id]
            time.sleep(sti_dur)
            yield "PuffOff", directions[trial_id]

            yield 'RegisterSti'
            yield 'RegisterBehavior'
        time.sleep(self.final_delay)
        yield 'RegisterBehavior'


# Odd Ball Test
class OBT:
    def __init__(self):
        print("Module [Odd Ball Test] loading...")
        self.initial_delay = 60
        print(f"Delay before first trial starts: {self.initial_delay} seconds")
        self.final_delay = 30
        print(f"Delay after last trial ends: {self.final_delay} seconds")
        self.redundant = 2
        self.deviant = 0
        print(f"Redundant Direction: {self.redundant}, Deviant Direction: {self.deviant}")
        self.inter_puff = 4
        print(f"Sessions are separated by {self.inter_puff} redundant puff")
        self.session_puff = 4
        print(f"Each session contains {self.session_puff} puff, one of them is deviant")
        self.inter_trial_delay = 10
        print(f"Duration between two trial varied between: {self.inter_trial_delay} [seconds]")

    def run(self):
        yield 'ShortPulse'
        time.sleep(self.initial_delay)

        t_list = []
        for _ in range(6):
            t_list += [self.redundant for _ in range(self.inter_puff)]
            t_list += list(np.random.permutation([self.deviant, ] +
                                                 [self.redundant for _ in range(self.session_puff-1)]))
        print(t_list)
        n_trial = len(t_list)
        for trial_id in range(n_trial):
            yield 'CheckCamera'
            puff_dir = t_list[trial_id]
            vis(f"{trial_id+1}/{n_trial}, Delay {self.inter_trial_delay}s, Direction: {puff_dir}")

            time.sleep(self.inter_trial_delay)
            if puff_dir == self.deviant:
                vis(f"Deviant Puff {self.deviant}", "G")
            else:
                vis(f"Redundant Puff {self.redundant}", "B")

            yield "PuffOn", puff_dir
            time.sleep(sti_dur)
            yield "PuffOff", puff_dir

            yield 'RegisterSti'
            yield 'RegisterBehavior'
        time.sleep(self.final_delay)
        yield 'RegisterBehavior'


# Odd Ball Control
class OBC:
    def __init__(self):
        print("Module [Odd Ball Control] loading...")
        self.initial_delay = 60
        print(f"Delay before first trial starts: {self.initial_delay} seconds")
        self.final_delay = 30
        print(f"Delay after last trial ends: {self.final_delay} seconds")
        self.redundant = 2
        self.deviant = 0
        print(f"Redundant Direction: {self.redundant}, Deviant Direction: {self.deviant}")
        self.inter_trial_delay = 10
        print(f"Duration between two trial varied between: {self.inter_trial_delay} [seconds]")

    def run(self):
        yield 'ShortPulse'
        time.sleep(self.initial_delay)
        #  0, 2, 4, 6
        directions = list(np.concatenate([np.random.permutation(8) for _ in range(6)]))
        print(directions)

        n_trial = len(directions)
        for trial_id in range(n_trial):
            yield 'CheckCamera'
            puff_dir = directions[trial_id]
            vis(f"{trial_id+1}/{n_trial}, Delay {self.inter_trial_delay}s, Direction: {puff_dir}")

            time.sleep(self.inter_trial_delay)
            if puff_dir == self.deviant:
                vis(f"Deviant Puff {self.deviant}", "G")
            elif puff_dir == self.redundant:
                vis(f"Redundant Puff {self.redundant}", "B")
            else:
                vis(f"Control Puff {puff_dir}", "R")
            yield "PuffOn", puff_dir
            time.sleep(sti_dur)
            yield "PuffOff", puff_dir

            yield 'RegisterSti'
            yield 'RegisterBehavior'
        time.sleep(self.final_delay)
        yield 'RegisterBehavior'


if __name__ == "__main__":
    print("-"*8+"Testing Modules"+"-"*8)
    x = OBC()
    start_time = time.time()
    for answer in x.run():
        print(time.time()-start_time, answer)

