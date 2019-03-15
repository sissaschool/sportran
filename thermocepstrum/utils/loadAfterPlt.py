

#ugly trick to posticipate the loading of the library
#import matplotlib.pyplot as plt


def loadRedefineGlobalPlt():
    global plt
    import matplotlib.pyplot as plt2
    plt=plt2

class Plt():
    def plot(self,*args, **kwargs):
        loadRedefineGlobalPlt()
        return plt2.plot(*args,**kwargs)

    def subplots(self,*args, **kwargs):
        loadRedefineGlobalPlt()
        return plt2.subplots(*args,**kwargs)

    def subplots_adjust(self,*args, **kwargs):
        loadRedefineGlobalPlt()
        return plt2.subplots_adjust(*args,**kwargs)


try:
   plt
except:
   plt=Plt()



