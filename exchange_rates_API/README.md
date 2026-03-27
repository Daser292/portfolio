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

2. Запустите сервис:

python app/main.py
Запуск через Docker

docker-compose up -d

### Запрос в БД

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

### Логирование

Ошибки записываются в файл logs/errors.log

### Названия переменных

DB_HOST	- Хост PostgreSQL
DB_PORT	- Порт PostgreSQL
DB_NAME	- Имя базы данных
DB_USER	- Пользователь БД
DB_PASSWORD	- Пароль БД
API_KEY	- Ключ API
API_URL	- URL API
FETCH_INTERVAL - Интервал запросов (минуты)
