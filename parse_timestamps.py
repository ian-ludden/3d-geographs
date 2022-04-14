from datetime import datetime, timedelta

if __name__ == '__main__':
    with open('./results/timestamps.txt', 'r') as f:
        dt_prev = None
        for line in f:
            dt_curr = datetime.strptime(line.strip(), "%m/%d/%y %H:%M:%S")
            if dt_prev:
                print((dt_curr - dt_prev).seconds + (dt_curr - dt_prev).days * 60 * 60 * 24)
            dt_prev = dt_curr
