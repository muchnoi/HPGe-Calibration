#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Я завел следующие поля в базе:
/EMS/DT /EMS/E /EMS/DE /EMS/B /EMS/DB /EMS/S /EMS/DS
Я добавил также простой "replace" -- работать будет только если время совпадает в точности.  
Время передается либо в time.time_struct либо float либо int (последние 2 -- секунды от 1970 г.).
"""
import sys, time, MySQLdb

def message(m):
    sys.stderr.write(str(m)+"\n")

def normalize(name):
    list = [ p for p in name.split("/") if p!="" ]    # remove spare '/'
    return "/"+"/".join(list)

def timeToMks(t):
    if isinstance(t, float) or isinstance(t, int) or isinstance(t, long):  return int(t*1000000)
    elif isinstance(t, time.struct_time):                                  return int(mktime(t)*1000000)
    else:                                                                  raise AttributeError("wrong type of time: (%s) %s"%(type(t),t))

queryName = "select c.c_id, c.c_type from t_entry e, t_domain d, t_channel c where e.c_channel=c.c_id and e.c_parent=d.c_id and d.c_path=%s and e.c_name=%s"

class Registry(object):

    def __init__(self, storage):
        self.storage  = storage;        self.registry = { }

    def registerChannel(self, name):
        if name in self.registry: raise KeyError(name)
        channel = Channel(name, self.storage.connection)
        if channel.id == None:    raise KeyError(name)  # could not obtain channel by name
        self.registry[name] = channel
        return channel

    def __call__(self, name):      return self.registry[name]
    def __contains__(self, name):  return name in self.registry
    def __iter__(self):            return iter(self.registry)

    def __getitem__(self, name):
        if not name in self.registry: raise KeyError(name)
        return self.registry[name].datum

    def __delitem__(self, name):
        if not name in self.registry: raise KeyError(name)
        del self.registry[name]

    def __setitem__(self, name, value):
        if not name in self.registry: raise KeyError(name)
        self.registry[name].modify(value)

class Channel(object):
    def __init__(self, name, connection):
        list = [ p for p in name.split("/") if p!="" ] # remove spare '/'
        self.id    = None
        self.stype = None
        self.name  = list[-1]
        self.path  = '/'+"/".join(list[:-1])
        assert self.name != ""
        assert self.path != ""
        self.datum    = None
        self.modified = False
        cursor = connection.cursor()
        cursor.execute(queryName, (self.path, self.name) )
        # message("name cursor messages: %s" % (cursor.messages,) )
        result = cursor.fetchall()
        cursor.close()

        # message("name result: %s" % (result,))
        if result!=None and len(result)>0:
            assert len(result)==1
            self.id    = result[0][0]
            self.stype = result[0][1]
            # message("%s registered as %d"%(name, self.id))

    def modify(self, value):
        assert self.id!=None
        assert self.stype in (0, 1)
        if self.stype==0: self.datum =   int(value)
        else:             self.datum = float(value)
        self.modified = True

class Storage(object):
    def __init__(self,\
                 host =   "sndfs1.sndonline",
                 user =   "author",
                 passwd = "power",
                 db     = "temdbase",
                 port   = 3307):
        self.connection = MySQLdb.connect(host=host, user=user, passwd=passwd, db=db, port=port)
        self.connection.autocommit(True)
        # message("connection messages: %s"%self.connection.messages)
        self.registry     = Registry(self)
        self.newstamp     = None
        self.oldstamp     = None

    def ping(self, reconnectFlag = False):    return self.connection.ping(reconnectFlag)

    def channel(self, name):                  return self.registry(name)

    def prepareQueues(self):
        queueD, queueI = [], []
        for name in self:
            channel = self.channel(name)
            if not channel.modified: continue
            if channel.stype==0: queueI.append(channel)
            else:                queueD.append(channel)
        return queueD, queueI

    def replace(self, t):
        self.newstamp = timeToMks(t)
        queueD, queueI = self.prepareQueues()
        if len(queueD)==0 + len(queueI)==0:
            message("nothing to replace");  return
        self.queries = self.generateReplace(queueD, "double") + self.generateReplace(queueI, "integer")
        self.execute()

    def insert(self, t):
        self.newstamp = timeToMks(t)
        queueD, queueI = self.prepareQueues()
        if len(queueD)==0 + len(queueI)==0:
            message("nothing to save");    return
        self.queries = [ self.generateInsert(queueD, "double"), self.generateInsert(queueI, "integer"), ]
        self.execute()

    def execute(self):
        self.connection.autocommit(False)
        cursor = self.connection.cursor()
        for q in self.queries:
            if q != None: cursor.execute(q)
#              message("execute: %s"%q)
        self.connection.commit()
        self.connection.autocommit(True)
        cursor.close()
        self.oldstamp = self.newstamp
        for name in self:
            channel = self.channel(name)
            if channel.modified: channel.modified = False

    def generateInsert(self, queue, table):
        if len(queue)<1: return None
        q     = "INSERT INTO t_data%s VALUES"%table
        space = " "
        for channel in queue:
            q += space
            q += "('%d', '%d', '%s')"%(channel.id, self.newstamp, channel.datum)
            space = ", "
        q += ";"
        return q

    def generateReplace(self, queue, table):
        result = [ ]
        if len(queue)<1: return result
        q = "UPDATE t_data%s SET c_value='%%s' where c_channel='%%d' and c_stamp='%%d';"%table
        for channel in queue:  result.append( q % (channel.datum, channel.id, self.newstamp) )
        return result

    def register(self, name0):
        name = normalize(name0); channel = self.registry.registerChannel(name);  return channel

    def __delitem__(self, name):        self.unregister(name, True)

    def unregister(self, name0, silent = False):
        name = normalize(name0);        del self.registry[name]

    def __getitem__(self, name0):
        name = normalize(name0);        return self.registry[name]

    def __setitem__(self, name0, value):
        name = normalize(name0);        self.registry[name] = value

    def __contains__(self, name0):
        name = normalize(name0);        return name in self.registry

    def __iter__(self):        return iter(self.registry)

if __name__ == "__main__" :
    import time

    # create control object (use temdbRO.Storage when importing this package !!)
    storage = Storage()

    # register ll channels of interesta
    dtime = storage.register("/TESTS/KOROL/DTIME")   # return channel -- for optimization
    storage.register("/TESTS/KOROL/ENERGY")          # plain registration
    storage.register("/TESTS/KOROL/DENERGY")         # plain registration

    for attempt in xrange(5):
      # put actual values/stamps
      dtime.modify(attempt) # same as next commented line, but dont use "dtime = attempt" !!:
#     storage["/TESTS/KOROL/DTIME"] = attempt

      energy = float((attempt+7)*999000 % 3450)/10.0
      storage["/TESTS/KOROL/ENERGY"]  = energy
      from math import sqrt
      storage["/TESTS/KOROL/DENERGY"] = sqrt( energy )

      tstamp = time.time()
      storage.insert( tstamp )

      dtime.modify((attempt+3)*10)
      storage["/TESTS/KOROL/DENERGY"] = 3*sqrt( energy )
      storage.replace( tstamp )

      if attempt!=4: time.sleep(1.3)
