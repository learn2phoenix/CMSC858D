import math
import numpy as np


iso_lens = [[100, 0, 50], [100, 100, 50]]
theta = [.5, .5]
points = [159, 28, 68]


def calc_likelihood():
    likelihood = 0
    for i in range(len(iso_lens)):
        iso_frac = theta[i] / sum(iso_lens[i])
        for j, exon in enumerate(iso_lens[i]):
            if exon == 0 or points[j] == 0:
                continue
            n_points = points[j]
            point_prob = 1 / n_points
            for exon_pos in range(exon):
                for _ in range(n_points):
                    likelihood += math.log(iso_frac * point_prob)
    return likelihood


def calc_denom():
    denom = []
    for exon_no, point_cnt in enumerate(points):
        for _ in range(point_cnt):
            point_prob = 0
            for i in range(len(iso_lens)):
                iso_frac = theta[i] / sum(iso_lens[i])
                for j, exon in enumerate(iso_lens[i]):
                    if exon == 0 or j != exon_no:
                        continue
                    for exon_pos in range(exon):
                        point_prob += iso_frac * (1 / point_cnt)
            denom.append(point_prob)
    return denom


def e_step():
    denom = calc_denom()
    expectation = [[], []]
    for i in range(len(iso_lens)):
        iso_frac = theta[i] / sum(iso_lens[i])
        for j, exon in enumerate(iso_lens[i]):
            for exon_no, point_cnt in enumerate(points):
                if exon == 0 or j != exon_no:
                    continue
                for point in range(point_cnt):
                    if exon_no > 0:
                        point += sum(points[:exon_no])
                    point_prob = 0
                    for exon_pos in range(exon):
                        point_prob += iso_frac * ((1 / point_cnt) / denom[point])
                    expectation[i].append(point_prob)
    return expectation


def m_step():
    expectation = e_step()
    return [sum(expectation[0]) / sum(points), sum(expectation[1]) / sum(points)]


step = 0
while step < 1000:
    calc_likelihood()
    theta = m_step()
    # print(f'Theta: {theta}')
    if step == 0:
        prev_likelihood = calc_likelihood()
    if step > 0:
        curr_likelihood = calc_likelihood()
        if abs(prev_likelihood - curr_likelihood) < 1e-3:
            break
        prev_likelihood = curr_likelihood
    step += 1

print(f'Finished in {step} steps. Theta: {theta}')


