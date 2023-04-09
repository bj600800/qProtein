"""
# ------------------------------------------------------------------------------
# Author:    Zhixin Dou
# Email:     bj600800@gmail.com
# DATE:      2023/03/06

# Description: Mongodb python interface on HPC server
# ------------------------------------------------------------------------------
"""


import os
import sys
import pymongo
from bson.objectid import ObjectId
import json
from gridfs import GridFS
import bson.binary


def printJson(object):
    jsonf = json.dumps(object, sort_keys=True, indent=4, separators=(',', ':'))
    print(jsonf)


def isChinese(ustring):
    for item in ustring:
        try:
            if u'\u4e00' <= item <= u'\u9fff':
                return True
        except UnicodeDecodeError:
            return True
    return False


class JSONEncoder(json.JSONEncoder):
    def default(self, o):
        if isinstance(o, ObjectId):
            return str(o)
        return json.JSONEncoder.default(self, o)


def objectIDToUnicode(id):
    if isinstance(id, ObjectId):
        #str = JSONEncoder().encode(id)
        return id.decode('utf8')
        # return str(id)
    else:
        return id


class MolService():
    def __init__(self):
        self.client = pymongo.MongoClient(host='211.87.224.168', port=4200,
                                          username='rcdb10', password='rcdb1234', authSource='admin')
        self.db = self.client['Mol-DB']
        self.fileDB = self.client['Mol-DB']
        # self.client.admin.authenticate('rcdb10', 'rcdb1234', mechanism='SCRAM-SHA-1', source='admin')
        self.user_Collection = self.db.user
        self.protA_Collection = self.db.A_prot_data
        self.protB_Collection = self.db.B_prot_data
        self.yeast = self.db.yeast_data_copy
        self.fs = GridFS(self.db, collection='pdbFiles')
        # self.cazy_gh = self.db.cazy_gh_alpha_beta_8
        # self.tibet = self.db.tibet

    def register(self, userName, password):
        user_Collection = self.db.user
        user_Collection.insert({'userName': userName, 'password': password})

    def login(self, userName, password):
        print('MolService--login')
        result = self.user_Collection.find_one({'userName': userName})
        if result and result['password'] == password:
            return True
        else:
            return False

    def putSeq(self, fasta_name, fasta_seq, mongo_collect):
        """
        update: 2022-10-19
        """
        str_colletc = 'self.db.'+mongo_collect
        collect = eval(str_colletc)

        condition = {'fasta_name': fasta_name}
        result = collect.find_one(condition)
        if result:
            result['fasta_name'] = fasta_name
            result['fasta_seq'] = fasta_seq
            record = collect.update_one(condition, {'$set': result})
        else:
            record = collect.insert_one(
                {'fasta_name': fasta_name, 'fasta_seq': fasta_seq})
        return record

    def getAllData(self, kwargs_dict, mongo_collect):
        str_colletc = 'self.db.' + mongo_collect
        collect = eval(str_colletc)

        # batch_size helps in case of the server_annot_API processing out of 10 min.
        return collect.find(kwargs_dict, batch_size=5)

    def setItem(self, id, key, value, mongo_collect):
        str_colletc = 'self.db.' + mongo_collect
        collect = eval(str_colletc)
        return collect.update_one({'_id': id}, {'$set': {key: value}})

    def setMultiItem(self, id, update_dict, mongo_collect):
        str_colletc = 'self.db.' + mongo_collect
        collect = eval(str_colletc)
        return collect.update({'_id': id}, {'$set': update_dict})

    def getProtA(self, seqName):
        return self.protA_Collection.find({'seqList.seqName': {'$eq': seqName}}, {'seqList.seqName.$': 1})

    def getProtB(self, seqName):
        return self.protB_Collection.find({'seqList.seqName': {'$eq': seqName}}, {'seqList.seqName.$': 1})

    def getFileInfo(self, fileID):
        result = self.fs.find_one({"_id": ObjectId(fileID)})
        return result

    def getFile(self, localFile, fileID):
        file = self.fs.find_one({"_id": ObjectId(fileID)})
        total = file.length
        inc = 0
        with open(localFile, "wb") as f:
            with self.fs.get(ObjectId(fileID)) as fd:
                for line in fd:
                    f.write(line)
                    inc = inc + len(line)
                    progress = int(float(inc) / float(total) * 100.0)
        f.close()
        return True

    def getBuffer(self, fileID):
        return self.fs.get(ObjectId(fileID)).read()

    def putFile(self, uploadFile):
        filename = os.path.basename(uploadFile)
        total = os.path.getsize(uploadFile)
        inc = 0
        with open(uploadFile, "rb") as f:
            buffer = f.read()
            invoice = self.fs.put(buffer, encoding='utf-8', filename=filename)
            inc = inc + len(buffer)
            progress = int(float(inc) / float(total) * 100.0)
        return invoice
