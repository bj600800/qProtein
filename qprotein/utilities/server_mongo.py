#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    mongo.py
# @Author:      Eric Dou
# @Time:        2022/10/11 17:00

"""
Mongodb python interface on HPC server

update: 2022-10-11
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

    def getBuffer(self, fileID):
        return self.fs.get(ObjectId(fileID)).read()

    def putFile(self, fileContent, filename, file_collect):
        str_colletc = 'self.db.' + file_collect
        collect = eval(str_colletc)

        condition = {'file_name': filename + '.cif'}
        result = collect.find_one(condition)
        if result:
            result['file_name'] = filename + '.cif'
            result['pdb_file'] = fileContent
            ret = collect.update_one(condition, {'$set': result})
            objid = result['_id']
        else:
            ret = collect.insert_one(
                {'file_name': filename + '.cif', 'pdb_file': fileContent})
            objid = ret.inserted_id

        return str(objid)
