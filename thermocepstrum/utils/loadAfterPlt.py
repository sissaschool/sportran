# ugly trick to posticipate the loading of the library
# import matplotlib.pyplot as plt


def loadRedefineGlobalPlt():
    global plt
    import matplotlib.pyplot as plt2_
    plt = plt2_


class Plt:

    def plot(self, *args, **kwargs):
        loadRedefineGlobalPlt()
        return plt.plot(*args, **kwargs)

    def subplots(self, *args, **kwargs):
        loadRedefineGlobalPlt()
        return plt.subplots(*args, **kwargs)

    def subplots_adjust(self, *args, **kwargs):
        loadRedefineGlobalPlt()
        return plt.subplots_adjust(*args, **kwargs)


try:
    plt
except:
    plt = Plt()
