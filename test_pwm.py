import sys
sys.path.insert(0, "/Users/Joshua Baller/Documents/Seelig/protease-sites/")
import compare_peptides as cp
sd = cp.scoring_distance()

pwm1 = cp.PWM('ACDR',sd)
(pwm1+pwm1)._curr_PWM

pwm2 = cp.PWM('CDRA',sd)
(pwm1+pwm2)._curr_PWM
(pwm1+pwm2+pwm2)._curr_PWM

pwm3 = cp.PWM('TTAC',sd)
(pwm1+pwm2+pwm2+pwm3)._curr_PWM