import sqlite3

class table:
    def __init__(self, db, cursor, table):
        self.db = db
        self.cursor = cursor
        self.table = table
        self.create_table(table)

    """ creates a database table if it doesn't already exist """
    def create_table(self, table):
        self.cursor.execute('''
            CREATE TABLE IF NOT EXISTS {} (id INTEGER PRIMARY KEY, str TEXT, mul INTEGER)
            '''.format(self.table))
        self.db.commit()

    """ creates a record, or increases multiplicity if it exists already """
    def add(self, id, str):
        if self.has(str):
            self.cursor.execute('''
                UPDATE {}
                SET mul = mul + 1
                WHERE str = ?
                '''.format(self.table),[str])
        else:
            self.cursor.execute('''
                INSERT INTO {} VALUES (?,?,?)
                '''.format(self.table), [id, str, 1])

    """ return a kmer based on its id """
    def get_str(self, id):
        return self.cursor.execute('''SELECT * FROM {} WHERE id = ?'''.format(self.table),
            [id]).fetchone()[1]
    """ return an id based on the kmer """
    def get_id(self, str):
        return self.cursor.execute('''SELECT * FROM {} WHERE str = ?'''.format(self.table),
            [str]).fetchone()[0]

    """ return whether or not a record exists in the database """
    def has(self, str):
        return self.cursor.execute('''SELECT EXISTS
            (SELECT 1 FROM {} WHERE str = ?)
            '''.format(self.table),[str]).fetchone()[0]

    """ print all the database records for debugging purposes """
    def print_all(self):
        recs = self.cursor.execute('''SELECT * FROM {}'''.format(self.table))
        for rec in recs:
            print (rec)


class db:
    def __init__(self, path):
        self.path = path

        self.db = sqlite3.connect(path)
        self.cursor = self.db.cursor()

        """ need two separate databases for nodes and edges for a graph """
        self.nodes = table(self.db, self.cursor, "nodes")
        self.edges = table(self.db, self.cursor, "edges")
