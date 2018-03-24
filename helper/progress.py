import time


def progress_bar(label, top, width):
    for i in range(1, top + 1):
        prog = round(i/top * width)
        print(("\r{label}: [{prog: <"+str(width)+"}] {fin}/{total}").format(
            label=label,
            prog=("#"*prog),
            fin=i,
            total=top
        ), end="")
        time.sleep(0.1)
    print("")


progress_bar("test", 100, 20)
