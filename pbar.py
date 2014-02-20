from __future__ import print_function
import sys, time

class ProgressBar:
    """
    Animated progresbar. Should work in console as well as the IPython
    notebook.
    """
    def __init__(self, its):
        self.its = its
        self.bar = '[]'
        self.sign = ':'
        self.width = 80
        self.amount(0)

    def update(self, it):
        print('\r', self.bar, end='')
        sys.stdout.flush()
        self.iteration(it + 1)

    def iteration(self, elapsed_it):
        self.amount((elapsed_it / float(self.its)) * 100.0)
        self.bar += '  %d of %s complete' % (elapsed_it, self.its)

    def amount(self, new_amount):
        pct_done = int(round((new_amount / 100.0) * 100.0))
        full_width = self.width - 2
        num_signs = int(round((pct_done / 100.0) * full_width))
        self.bar = '[' + self.sign * num_signs + ' ' * (full_width - num_signs) + ']'
        pct_place = (len(self.bar) // 2) - len(str(pct_done))
        pct_string = '%d%%' % pct_done
        self.bar = self.bar[0:pct_place] + \
            (pct_string + self.bar[pct_place + len(pct_string):])
