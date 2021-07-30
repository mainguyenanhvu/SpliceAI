#logging_config.yml

version: 1
disable_existing_loggers: False
formatters:
  standard:
    format: '%(asctime)s [%(levelname)s] %(name)s - %(message)s'
handlers:
  console:
    class: logging.StreamHandler
    level: INFO
    formatter: standard
    stream: ext://sys.stdout
  file:
    class: logging.FileHandler
    level: WARNING
    formatter: standard
    filename: spliceai.log
  email:
    class: logging.handlers.SMTPHandler
    level: WARNING
    mailhost: smtp.gmail.com
    fromaddr: to@address.co.uk
    toaddrs: vumna@vinbigdata.org
    subject: Oh no, something's gone wrong!
    credentials: [email, password]
    secure: []
loggers:
  cloaked_chatter:
    level: INFO
    handlers: [console, file]
    propagate: no
  root:
    level: WARNING
    handlers: [console, file, email]
    propagate: True