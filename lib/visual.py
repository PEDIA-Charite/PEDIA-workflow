import multiprocessing
from functools import wraps, partial


def print_status(label, width, cur, size):
    '''Print a statusbar.'''
    len_size = str(len(str(size)))
    # prevent division by zero by defaulting to 100% progress
    prog = round(cur/size * width) if size > 0 else 1
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
    def progress_decorator(mapped_func):
        @wraps(mapped_func)
        def progress_wrapper(iterable, *args, **kwds):
            result = []
            size = len(iterable)

            for i, item in enumerate(iterable):
                print_status(label, width, i+1, size)
                result.append(mapped_func(item, *args, **kwds))
            print("")
            return result
        return progress_wrapper

    return progress_decorator


def multiprocess(
        label, map_func, iterable, *args, **kwargs
):
    with multiprocessing.Pool() as pool:
        result = []
        size = len(iterable)
        for i, res in enumerate(
                pool.imap_unordered(
                    partial(map_func, *args, **kwargs), iterable
                )
        ):
            print_status(label, 20, i+1, size)
            result.append(res)
        print("")
    return result


if __name__ == '__main__':
    import time

    @progress_bar("Test")
    def test(iterable, test):
        for i in iterable:
            time.sleep(1)
            yield i + test
    test(range(10), 100)
