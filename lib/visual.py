from functools import wraps


def print_status(label, width, cur, size):
    '''Print a statusbar.'''
    len_size = str(len(str(size)))
    prog = round(cur/size * width)
    print(
        ("\r{label}: \t[{prog: <"
         + str(width)
         + "}] {fin:>"+len_size+"}/{total:>"+len_size+"}").format(
             label=label,
             prog=("#"*prog),
             fin=cur,
             total=size
         ),
        end="")


def progress_bar(label, width=20):
    '''Build progressbar around yielding function.'''
    def progress_decorator(iter_func):
        @wraps(iter_func)
        def progress_wrapper(*args, **kwds):
            result = []
            size = len(args[0])
            for i, res in enumerate(iter_func(*args, **kwds)):
                print_status(label, width, i+1, size)
                result.append(res)
            print("")
            return result
        return progress_wrapper

    return progress_decorator


if __name__ == '__main__':
    import time

    @progress_bar("Test")
    def test(iterable, test):
        for i in iterable:
            time.sleep(1)
            yield i + test
    test(range(10), 100)
