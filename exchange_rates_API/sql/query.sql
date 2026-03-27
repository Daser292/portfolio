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