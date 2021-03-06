#from progressbar import *
#from blessings   import Terminal

#def pbar(title, *loop_length):

#    if loop_length:     pb = ProgressBar(widgets = pbar_widgets(title), \
#                        term_width = terminal_width(), \
#                        maxval = loop_length)

#    if not loop_length: pb = ProgressBar(widgets = pbar_widgets(title), \
#                        term_width = terminal_width())

#    return pb

#def pbar_widgets(title): return [title + ': ',                               \
#                                 Percentage(),                               \
#                                 ' ',                                        \
#                                 Bar(marker = '-', left = '|', right = '|'), \
#                                 ' ', ETA(), ' ']

#def term_width(): return Terminal().width

def term_width():

    import os

    env = os.environ

    def ioctl_GWINSZ(fd):
        try:
            import fcntl, termios, struct, os
            cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,
        '1234'))
        except:
            return
        return cr
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        cr = (env.get('LINES', 25), env.get('COLUMNS', 80))

        ### Use get(key[, default]) instead of a try/catch
        #try:
        #    cr = (env['LINES'], env['COLUMNS'])
        #except:
        #    cr = (25, 80)
#    return int(cr[1]), int(cr[0])
    return int(cr[1])
