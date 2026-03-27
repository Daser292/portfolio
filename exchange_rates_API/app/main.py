import os
import time
import schedule
from dotenv import load_dotenv
from database import Database
from fetcher import CurrencyFetcher
from logger_config import setup_logger

# Загружаем настройки из файла .env
load_dotenv()

# Настраиваем логирование
info_logger, error_logger = setup_logger()

def fetch_and_save():
    """Главная функция: получить курсы и сохранить в базу"""
    
    info_logger.info("=" * 50)
    info_logger.info("Начинаю получение курсов валют")
    
    # Создаем объекты для работы с базой и API
    db = Database()
    fetcher = CurrencyFetcher(
        api_key=os.getenv('API_KEY'),
        api_url=os.getenv('API_URL')
    )
    
    # Подключаемся к базе данных
    if not db.connect():
        error_logger.error("Не удалось подключиться к базе данных")
        return
    
    # Создаем таблицы, если их еще нет
    if not db.create_tables():
        error_logger.error("Не удалось создать таблицы")
        db.close()
        return
    
    # Получаем курсы валют
    result = fetcher.fetch_rates()
    
    if result['success']:
        info_logger.info(f"Успешно получены данные. Время ответа: {result['response_time']}мс")
        
        # Сохраняем информацию о запросе
        request_id = db.insert_request(
            status='success',
            response_time_ms=result['response_time']
        )
        
        if request_id:
            # Извлекаем курсы из ответа API
            rates = fetcher.parse_rates(result['data'])
            
            # Сохраняем каждый курс в базу
            for rate in rates:
                if db.insert_response(
                    request_id=request_id,
                    currency_from=rate['from'],
                    currency_to=rate['to'],
                    rate=rate['rate']
                ):
                    info_logger.info(f"Сохранен курс: {rate['from']} -> {rate['to']} = {rate['rate']}")
                else:
                    error_logger.error(f"Не удалось сохранить курс: {rate['from']} -> {rate['to']}")
        else:
            error_logger.error("Не удалось сохранить запись о запросе")
            
    else:
        # Если произошла ошибка при получении данных
        error_logger.error(f"Ошибка при получении данных: {result['error']}")
        
        # Сохраняем информацию об ошибке
        db.insert_request(
            status='error',
            response_time_ms=result['response_time'],
            error_message=result['error']
        )
    
    # Закрываем соединение с базой
    db.close()
    info_logger.info("Завершение получения курсов")
    info_logger.info("=" * 50)

def main():
    """Главная функция, запускает программу"""
    
    # Получаем интервал из настроек (по умолчанию 10 минут)
    interval = int(os.getenv('FETCH_INTERVAL', 10))
    
    info_logger.info(f"Сервис запущен. Интервал: {interval} минут")
    
    # Запускаем сразу при старте
    fetch_and_save()
    
    # Настраиваем расписание
    schedule.every(interval).minutes.do(fetch_and_save)
    
    # Бесконечный цикл, который ждет следующего запуска
    while True:
        schedule.run_pending()
        time.sleep(1)

# Точка входа в программу
if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        info_logger.info("Сервис остановлен пользователем")
    except Exception as e:
        error_logger.error(f"Критическая ошибка: {str(e)}")
        raise