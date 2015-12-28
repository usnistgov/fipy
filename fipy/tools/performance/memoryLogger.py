#!/usr/bin/env python

import os
import select
import signal
import threading
import time

__all__ = ["MemoryHighWaterThread", "MemoryLogger"]

class MemoryHighWaterThread(threading.Thread):
    def __init__(self, pid, sampleTime = 1):
        threading.Thread.__init__(self)
        self.pid = pid
        self.sampleTime = sampleTime
        self.actionDone = threading.Event()
        self.setDaemon(True)

        self.maxMem = 0

    def run(self):
        self.maxMem = 0

        while not self.actionDone.isSet() or self.maxMem == 0:
            ff = os.popen('ps -p %i -o vsz' % self.pid)
            ff.readline() # skip a line
            s = ff.readline()

            ff.close()

            self.maxMem = max(self.maxMem, int(s))
            time.sleep(self.sampleTime)

    def stop(self):
        self.actionDone.set()
        self.join()

        return self.maxMem

class MemoryLogger:
    def __init__(self, sampleTime = 1):
        self.pid = os.getpid()

        self.r, self.w = os.pipe()
        self.forkID = os.fork()

        if self.forkID == 0:
            os.close(self.r)

            self.thread = None

            def startHandler(signum, frame):
                if self.thread is not None:
                    self.stop()

                self.thread = MemoryHighWaterThread(pid=self.pid, sampleTime=sampleTime)
                self.thread.start()
##                 print "start"


            def stopHandler(signum, frame):
                if self.thread is not None:
                    maxMem = self.thread.stop()
                    self.thread = None
                else:
                    maxMem = -9999

                os.write(self.w, str(maxMem))

            def hupHandler(signum, frame):
                os._exit(1)

            signal.signal(signal.SIGUSR1, startHandler)
            signal.signal(signal.SIGUSR2, stopHandler)
            signal.signal(signal.SIGHUP, hupHandler)

            while 1:
                time.sleep(10000000)

        else:
            os.close(self.w)

    def __del__(self):
        import signal
        os.kill(self.forkID, signal.SIGHUP)

    def start(self):
        os.kill(self.forkID, signal.SIGUSR1)

    def stop(self):
        os.kill(self.forkID, signal.SIGUSR2)
        inputs, outputs, exceptions = select.select([self.r], [], [], 2)
        if self.r in inputs:
            maxMem = int(os.read(self.r, 1024))
        else:
            maxMem = -8888
        return maxMem


if __name__ == "__main__":
    print "MemoryHighWaterThread"
    for attempt in range(10):
        thread = MemoryHighWaterThread(pid=os.getpid(), sampleTime=1)
        thread.start()
        print thread.stop()

    print "MemoryLogger"
    logger = MemoryLogger(sampleTime=1)
    for attempt in range(10):
        logger.start()
        print logger.stop()
