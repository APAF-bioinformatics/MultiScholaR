"""KIT-262 M3 budget arithmetic (local, standard library only). Core-hours = job-minutes x 4 / 60, round(., 2); 10 job-minutes of set-up
per job. Inputs are the ticket's basis (section 11: 39.2 s Automatic journey, 17.8 s Manual journey at n = 12; stress cells twice that;
E10 7.2 s per stress unit; E8 2,850 s) until L3 replaces them with its own seconds (G2). Spread 1.93 = the largest ratio of a counted job's
sweep seconds to its lifecycle's seen in M2 (AUDIT-KIT-262-m1m2 Q3). Run: python3 -I docs/qa/workflow-execution/kit/kit-262/m3-budget-k262.py"""
ch = lambda m: round(m * 4 / 60, 2)
AUTO, MAN, E10 = 39.2, 17.8, 7.2
STRESS = 2 * (AUTO + MAN)  # P + C1 journeys at n = 24 (twice the n = 12 time)
SP = 1.93
USED = 87.24  # recomputed from runs.tsv by the audit (RC7): 58.56 through M2 and its record fixes + 28.68 re-measurement
CEIL = 250.32  # VAL-390
print("per-unit seconds: eight cells %.1f (P + C1), stress %.1f + E10 %.1f" % (AUTO + MAN, STRESS, E10))

# L3
eight_main = 16 * (AUTO + MAN + MAN)      # P, C1 Manual, forced E11 Manual
eight_e8 = 8 * (AUTO + MAN + MAN)
stress_main = 4 * (STRESS + E10 + MAN * 2)
stress_e8 = 2 * (STRESS + E10 + MAN * 2)
e9 = 16 * AUTO
s3 = eight_main + eight_e8 + stress_main + stress_e8 + e9
s5 = 10 * (STRESS + E10 + MAN * 2)
print("S3 CPU-s %.0f -> %.1f min at 4 workers (spread 1.0), %.1f (1.93)" % (s3, s3 / 240, s3 * SP / 240))
print("S5 CPU-s %.0f -> %.1f min, %.1f" % (s5, s5 / 240, s5 * SP / 240))
l3_exp = 1 + 5 + s3 / 240 + 2 + s5 / 240 + 10
l3_hi = 1 + 5 + s3 * SP / 240 + 3 + s5 * SP / 240 + 10
print("L3 expected %.1f job-min = %.2f core-h; at 1.93 %.1f = %.2f; ticket line 2.00 / 4.00" % (l3_exp, ch(l3_exp), l3_hi, ch(l3_hi)))
timeouts = 5 + 20 + 50 + 10 + 15
print("L3 timeouts sum %d script-min + 10 = %d job-min = %.2f core-h (limits for a hung script)" % (timeouts, timeouts + 10, ch(timeouts + 10)))

# stage 1: 100 stress units (S-N, S-E at 1:50), rows private
st1 = 100 * (STRESS + E10)
print("stage 1 %.0f CPU-s: expected %.1f job-min = %.2f; at 1.93 %.1f = %.2f" % (st1, st1 / 240 + 10, ch(st1 / 240 + 10), st1 * SP / 240 + 10, ch(st1 * SP / 240 + 10)))

# M3 counted (O-1 (ii), O-4 (a)); eight cells 1,550 units; stress at R = 250 + 200 (expected) or doubled (ceiling)
exp_cpu = 1550 * (AUTO + MAN) + 450 * (STRESS + E10) + 2850
ceil_cpu = SP * (1550 * (AUTO + MAN) + 900 * (STRESS + E10) + 2850)
exp_min = exp_cpu / 240 + 5 * 10
ceil_min = ceil_cpu / 240 + 12 * 10
print("M3 expected %.0f CPU-s -> %.2f job-min = %.2f core-h (ticket 657.25 / 43.82); ceiling at R doubled, 1.93: %.0f CPU-s -> %.1f = %.2f (ticket's widest 128.03)" %
      (exp_cpu, exp_min, ch(exp_min), ceil_cpu, ceil_min, ch(ceil_min)))
sc_left = round(2.67 - 0.80, 2)
exp_total = USED + ch(l3_exp) + ch(st1 / 240 + 10) + ch(exp_min) + 1.00
ceil_line = USED + 4.00 + ch(st1 * SP / 240 + 10) + ch(ceil_min) + sc_left
print("running line: used %.2f of %.2f, remaining %.2f" % (USED, CEIL, CEIL - USED))
print("expected total %.2f; after L3 at ceiling %.2f; after stage 1 at 1.93 %.2f; after M3 at ceiling %.2f; after scoring %.2f; with the reserves (8.77 + 7.43) %.2f; headroom %.2f" %
      (exp_total, USED + 4.00, USED + 4.00 + ch(st1 * SP / 240 + 10), USED + 4.00 + ch(st1 * SP / 240 + 10) + ch(ceil_min), ceil_line, ceil_line + 8.77 + 7.43,
       CEIL - (ceil_line + 8.77 + 7.43)))
# the planning figure of R above which the stress effect cell routes to the coordinator
print("R ~ 3,354 x SD^2: SD 0.265 -> %.0f, 0.35 -> %.0f, 0.40 -> %.0f; R > 400 above SD %.3f" % (3354 * 0.265 ** 2, 3354 * 0.35 ** 2, 3354 * 0.40 ** 2, (400 / 3354) ** 0.5))
