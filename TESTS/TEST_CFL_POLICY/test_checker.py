"""Negative controls for accepted-time/event evidence, not solver physics."""
import unittest
import numpy as np
from check_cases import check_times, DT0, MAX_DT, END, OUTPUT


class TimeCheckerTests(unittest.TestCase):
    @staticmethod
    def trajectory():
        time=0.0; old=older=DT0; event=OUTPUT; rows=[]
        while time < END:
            dt=min(MAX_DT,END-time,event-time,1.1*0.5*(old+older))
            rows.append([time,dt,MAX_DT])
            time+=dt; older,old=old,dt
            if abs(time-event)<1e-14:
                event=min(END,event+OUTPUT)
        return np.array(rows)

    def test_complete_trajectory(self):
        self.assertEqual(check_times(self.trajectory(),True)["output_hits"],3)

    def test_reject_bad_time_or_cap(self):
        for column,amount in ((0,1e-4),(1,1e-4),(2,MAX_DT)):
            rows=self.trajectory(); rows[2,column]+=amount
            with self.assertRaises(AssertionError): check_times(rows,True)

    def test_reject_nonfinite_or_missing_final(self):
        rows=self.trajectory(); rows[2,1]=np.nan
        with self.assertRaises(AssertionError): check_times(rows)
        with self.assertRaises(AssertionError): check_times(self.trajectory()[:-1])


if __name__ == "__main__":
    unittest.main()
