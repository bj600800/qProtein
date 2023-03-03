#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @Filename:    main.py
# @Author:      Eric Dou
# @Time:        2022/10/11 16:59

"""
The backends main script of pyMulstruct on HPC server.

Including:
    1.Summary results and upload to database.
    2.Communication between server and database.

Update: 2022-11-3
"""

import server_annot_API
import server_annot_Local
import server_mongo

# import server_mulstruct_tools


class ServerMulstruct:
    def __init__(self, service):
        self.service = service

    def update_annotation(self):
        """
        Update sequence annotation using annot_API, annot_third and mulstruct_tools.

        including: Protein function, protein structure from uniprot,
                    thermal stability, protein structure from local alphafold2,
                    CAZy family classification.
        """
        server_annot_API.AnnotCrawl(self.service)

    def main(self):
        self.update_annotation()


def run():
    """
    main function
    """
    service = server_mongo.MolService()
    obj = ServerMulstruct(service)
    obj.main()


run()
