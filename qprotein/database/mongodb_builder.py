import sqlite3
con = sqlite3.connect("DEMO.db")
cur = con.cursor()
sql = "CREATE TABLE IF NOT EXISTS test(id INTEGER PRIMARY KEY,name TEXT,age INTEGER)"
cur.execute(sql)
data = "1,'Desire',5"
cur.execute('INSERT INTO test VALUES (%s)' % data)
cur.execute("INSERT INTO test values(?,?,?)", (6, "zgq", 20))
cur.executemany('INSERT INTO test VALUES (?,?,?)', [
                (3, 'name3', 19), (4, 'name4', 26)])
cur.execute("select * from Test")
print(cur.fetchall())
