# NGS-emulator

Я думаю,что можно в этом проекте выделить три основных этапа.
На первом создать алгоритм и способ его фунционирования.
Вторым можно назвать этап разработки API приложения (интерфейса и способов взаимодействия с пользователем).
И третим этапом я бы выделил накручивание опционала, который у меня в голове уже вертиться.

Итак, сейчас нада состредоточиться на первом этапе - создание алгоритма, а по сути основного ядра программы.
Я вижу это как создание независимого работающего пакета в любом виде. Скорее всего, пакета Python, по понятным причинам.
Для декомпозиции предлагаю разбить эту программу на три блока а задачу эмуляции, соответственно на 3 подзадачи:

  1. **Создание первичного API**. По сути, в этом блоке, нужно определиться с деталями и границами
  нашего проекта косательно его функционала. Определить ЧТо мы хотим пользователю эмулировать,
  и какую информацию мы хотим узнать у пользователя для выполнения целей.
  
  2. **Создание ядра программ** а именно механизмов построения данных. Определиться как и в каком
  разрезе использовать рандом, какие распределения использовать для симуляции случайных величин,
  как ограничить случайность данных, для соблюдения требований пользователя к эмулируемому кейсу.
  
  3. **Создание конвертатора в файлы**. Полученные данные из второго блока пользователю необходимы
  в виде конкретных форматов (vcf, csv, tsv, bad, baf, logr и др.). Зелью данного блока можно
  назвать создание скрипта, который на основании запроса из 1-ого блока сформует данные из 2-ого
  блока и запишет фаил
  
Очевидно, что данные блоки несовсем независимы, однако с точки зрения, написания чернового варианта кода вполне
могут быть написанны по отдельности.

v0.1
Выложил лысый код, эмулирующий LogR и Baf файлы с количеством ридов (chr pos DP и AD поля) без заголовков.
От пользователя в качестве аргумента принимает размер данных (количество позиций) и количество
суб-клональных популяций.
  
