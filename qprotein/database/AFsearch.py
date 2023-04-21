"""
# ------------------------------------------------------------------------------
# Author:    Dou Zhixin
# Email:     bj600800@gmail.com
# DATE:      2023/04/21

# Description: 
# ------------------------------------------------------------------------------
"""
from qprotein.database.sqlite3_searcher import SqlSearch
from qprotein.utilities import logger

logger = logger.setup_log(name=__name__)


class AFmetaSearch(SqlSearch):
    def __init__(self, sql_db, table_name):
        super(AFmetaSearch, self).__init__(sql_db)

        self.table_name = table_name

    def search_sql(self, cursor, target_list):
        results = []
        for i in target_list:
            sql_query = f"SELECT * FROM {self.table_name} WHERE UniDesc=? Limit 1"
            cursor.execute(sql_query, (i[0],))
            _results = cursor.fetchall()
            results.append(_results)
        return results

    def write_output(self, results, file):
        with open(file, 'w') as wt:
            for result in results:
                for one_result in result:
                    str_result = [str(i) for i in one_result]
                    wt.write('\t'.join(str_result) + '\n')

    def run(self):
        cursor = self.connect_sql(self.sql_db)
        sql_target = f"SELECT UniDesc, COUNT (UniDesc) FROM {self.table_name} GROUP BY UniDesc ORDER BY COUNT (UniDesc) DESC LIMIT 26"
        logger.info(f"Fetch targets from {self.table_name}")
        target_list = self.fetch_results(cursor, sql_target)

        logger.info("Search results")
        results = self.search_sql(cursor, target_list)
        file = r'D:\subject\active\1-qProtein\data\alphafold_metadata\two-results_for_function_cluster.txt'

        logger.info(f"Write the results in {file}")
        self.write_output(results, file)
