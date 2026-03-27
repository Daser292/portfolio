import requests
import time

class CurrencyFetcher:
    def __init__(self, api_key, api_url):
        self.api_key = api_key
        self.api_url = api_url
    
    def fetch_rates(self):
        """Получение курсов валют с API"""
        start_time = time.time()
        
        try:
            # Отправляем запрос к API
            response = requests.get(
                self.api_url,
                params={'apikey': self.api_key},
                timeout=10  # ждем ответ максимум 10 секунд
            )
            
            # Считаем сколько миллисекунд занял запрос
            response_time = int((time.time() - start_time) * 1000)
            
            # Проверяем, успешно ли ответил сервер
            if response.status_code == 200:
                data = response.json()
                return {
                    'success': True,
                    'response_time': response_time,
                    'data': data
                }
            else:
                return {
                    'success': False,
                    'response_time': response_time,
                    'error': f"API вернул статус {response.status_code}"
                }
                
        except requests.exceptions.Timeout:
            response_time = int((time.time() - start_time) * 1000)
            return {
                'success': False,
                'response_time': response_time,
                'error': "Таймаут при запросе к API"
            }
        except requests.exceptions.ConnectionError:
            response_time = int((time.time() - start_time) * 1000)
            return {
                'success': False,
                'response_time': response_time,
                'error': "Ошибка подключения к API"
            }
        except Exception as e:
            response_time = int((time.time() - start_time) * 1000)
            return {
                'success': False,
                'response_time': response_time,
                'error': f"Неизвестная ошибка: {str(e)}"
            }
    
    def parse_rates(self, data):
        """Извлекает нужные курсы из ответа API"""
        rates = []
        base_currency = data.get('base', 'USD')
        
        # Берем курсы к рублю и евро для примера
        currencies = ['RUB', 'EUR']
        
        for currency in currencies:
            if currency in data.get('rates', {}):
                rates.append({
                    'from': base_currency,
                    'to': currency,
                    'rate': data['rates'][currency]
                })
        
        return rates