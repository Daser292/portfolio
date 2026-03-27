import logging
import os

def setup_logger():
    """Настройка логирования"""
    
    # Создаем директорию для логов если её нет
    if not os.path.exists('portfolio/exchange_rates_API/logs'):
        os.makedirs('portfolio/exchange_rates_API/logs')
    
    # Формат для логов
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Логгер для ошибок
    error_logger = logging.getLogger('error_logger')
    error_logger.setLevel(logging.ERROR)
    
    # Файловый обработчик для ошибок
    file_handler = logging.FileHandler(f'portfolio/exchange_rates_API/logs/errors.log', encoding='utf-8')
    file_handler.setLevel(logging.ERROR)
    file_handler.setFormatter(formatter)
    error_logger.addHandler(file_handler)
    
    # Логгер для обычной информации
    info_logger = logging.getLogger('info_logger')
    info_logger.setLevel(logging.INFO)
    
    # Консольный обработчик
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    info_logger.addHandler(console_handler)
    
    return info_logger, error_logger