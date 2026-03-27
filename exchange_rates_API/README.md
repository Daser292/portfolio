<<<<<<< HEAD
# Currency Fetcher Service

Сервис для получения курсов валют с публичного API и сохранения в PostgreSQL.

## Что делает этот сервис

1. Раз в N минут запрашивает курсы валют с API
2. Сохраняет результаты в базу данных PostgreSQL
3. Логирует все ошибки в отдельный файл
4. Позволяет посмотреть историю запросов через SQL

## Как запустить

### Локальный запуск

1. Создайте базу данных PostgreSQL:

createdb currency_db

Установите зависимости:

pip install -r requirements.txt

Скопируйте .env.example в .env и заполните настройки:

cp .env.example .env
Запустите сервис:

python app/main.py
Запуск через Docker

docker-compose up -d
SQL запрос для просмотра истории

sql
-- Выводит все запросы и полученные курсы
SELECT 
    r.id AS request_id,
    r.request_time,
    r.status,
    r.response_time_ms,
    r.error_message,
    res.currency_from,
    res.currency_to,
    res.rate,
    res.fetched_at
FROM requests r
LEFT JOIN responses res ON r.id = res.request_id
ORDER BY r.request_time DESC, res.fetched_at DESC;
Структура базы данных

requests — информация о запросах к API (время, статус, ошибки)
responses — полученные курсы валют (связаны с requests)
Логирование

Ошибки записываются в файл logs/errors.log

Переменные окружения

Переменная	Описание
DB_HOST	Хост PostgreSQL
DB_PORT	Порт PostgreSQL
DB_NAME	Имя базы данных
DB_USER	Пользователь БД
DB_PASSWORD	Пароль БД
API_KEY	Ключ API
API_URL	URL API
FETCH_INTERVAL	Интервал запросов (минуты)
text
=======
# Баулин Семен — Backend Python разработчик

## Обо мне

Backend разработчик с фокусом на разработку серверной логики, обработку данных и построение масштабируемых систем. Имею опыт работы с Python для создания API, обработки данных и автоматизации процессов.
Обладаю сильной базой в алгоритмах, работе с базами данных и написании чистого, поддерживаемого кода.

Владею полным стеком технологий для работы с данными:

- Языки: Python (Pandas, SciPy, Statsmodels), SQL (оконные функции, оптимизация запросов)
- Backend: FastAPI, Flask
- Визуализация: Plotly, Seaborn, Matplotlib, Tableau
- Базы данных: PostgreSQL, MySQL, Oracle
- Прочее: Git, Jupyter Notebook, Apache Airflow, Docker, Linux, REST API, написание ETL-процессов

## Проекты

### Diplom

Диплом на тему: Реализация метода сглаженных частиц на задачах гидродинамики

Разработал вычислительную часть на C++ для моделирования гидродинамических процессов
Реализовал интеграцию с Python для обработки результатов
Построил визуализацию временных итераций с использованием matplotlib
Организовал пайплайн обработки данных: генерация → сохранение → анализ → визуализация

- Технологии: Python, C++
- Код на C++ берет требуемую задачу гидродинамики, высчитывает результаты, записывает в папку /out
 далее через pyhon (matplotlib) выводятся результаты на каждой временной итерации

### Дополнительно

- Опыт работы с большими объёмами данных
- Оптимизация SQL-запросов и обработка данных
- Навыки автоматизации и построения воспроизводимых пайплайнов

## Контакты

- Email: vasisdas36@gmail.com
- Telegram: @sambaulin
- Резюме: https://hh.ru/resume/cefbd7a8ff0b5567340039ed1f34307845365a
>>>>>>> 5dd1a09b7c890a61d928653730b454a796aec496
