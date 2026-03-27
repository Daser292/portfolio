import psycopg2
import os
from dotenv import load_dotenv

load_dotenv(dotenv_path='portfolio/exchange_rates_API/.env')

class Database:
    def __init__(self):
        self.conn = None
        self.cursor = None
    
    def connect(self):
        """Подключение к базе данных"""
        try:
            self.conn = psycopg2.connect(
                host=os.getenv('DB_HOST'),
                port=os.getenv('DB_PORT'),
                database=os.getenv('DB_NAME'),
                user=os.getenv('DB_USER'),
                password=os.getenv('DB_PASSWORD')
            )
            self.cursor = self.conn.cursor()
            return True
        except Exception as e:
            print(f"Ошибка подключения: {e}")
            return False
    
    def create_tables(self):
        """Создание таблиц"""
        try:
            # Таблица requests - хранит информацию о каждом запросе
            self.cursor.execute("""
                CREATE TABLE IF NOT EXISTS requests (
                    id SERIAL PRIMARY KEY,
                    request_time TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    status VARCHAR(50),
                    response_time_ms INTEGER,
                    error_message TEXT
                )
            """)
            
            # Таблица responses - хранит полученные курсы валют
            self.cursor.execute("""
                CREATE TABLE IF NOT EXISTS responses (
                    id SERIAL PRIMARY KEY,
                    request_id INTEGER REFERENCES requests(id) ON DELETE CASCADE,
                    currency_from VARCHAR(10),
                    currency_to VARCHAR(10),
                    rate DECIMAL(10, 4),
                    fetched_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
                )
            """)
            
            self.conn.commit()
            return True
        except Exception as e:
            self.conn.rollback()
            print(f"Ошибка создания таблиц: {e}")
            return False
    
    def insert_request(self, status, response_time_ms, error_message=None):
        """Сохраняет запись о том, что был сделан запрос к API"""
        try:
            self.cursor.execute("""
                INSERT INTO requests (status, response_time_ms, error_message)
                VALUES (%s, %s, %s)
                RETURNING id
            """, (status, response_time_ms, error_message))
            self.conn.commit()
            return self.cursor.fetchone()[0]
        except Exception as e:
            self.conn.rollback()
            return None
    
    def insert_response(self, request_id, currency_from, currency_to, rate):
        """Сохраняет полученный курс валюты"""
        try:
            self.cursor.execute("""
                INSERT INTO responses (request_id, currency_from, currency_to, rate)
                VALUES (%s, %s, %s, %s)
            """, (request_id, currency_from, currency_to, rate))
            self.conn.commit()
            return True
        except Exception as e:
            self.conn.rollback()
            return False
    
    def close(self):
        """Закрывает соединение с базой данных"""
        if self.cursor:
            self.cursor.close()
        if self.conn:
            self.conn.close()