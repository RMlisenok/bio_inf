import paramiko
from paramiko import SSHClient, AutoAddPolicy
import time


class SSHTunnel:
    def __init__(self):
        self.clients = []

    def connect(self, host, port, username, password):
        client = SSHClient()
        client.set_missing_host_key_policy(AutoAddPolicy())

        if self.clients:
            # Используем последнее соединение как транспорт
            transport = self.clients[-1].get_transport()
            dest_addr = (host, port)
            local_addr = ('localhost', 22)
            channel = transport.open_channel("direct-tcpip", dest_addr, local_addr)
            client.connect(hostname=host, port=port, username=username,
                           password=password, sock=channel)
        else:
            # Первое соединение
            client.connect(host, port=port, username=username, password=password)

        self.clients.append(client)
        print(f"✓ Подключено к {username}@{host}:{port}")
        time.sleep(1)  # Пауза между подключениями
        return client

    def execute(self, commands):
        if not self.clients:
            raise Exception("Нет активных подключений")

        results = []
        for cmd in commands:
            stdin, stdout, stderr = self.clients[-1].exec_command(cmd)
            output = stdout.read().decode().strip()
            error = stderr.read().decode().strip()

            results.append({
                'command': cmd,
                'output': output,
                'error': error
            })
            print(f"\nCommand: {cmd}")
            print(f"Output:\n{output}")
            if error:
                print(f"Errors:\n{error}")

        return results

    def close(self):
        for client in reversed(self.clients):
            client.close()
        print("Все соединения закрыты")


# Конфигурация
JUMP_CHAIN = [
    ("37.252.1.93", 2003, "makashov", "UuhiWCd.fo-1sB{&+R]YA%)TwX^8~3"),
    ("127.0.0.1", 2217, "science2", "DL]M]0_HloCAJayh)oMpI$KQaXtl}w"),
    ("192.168.1.4", 22, "student1", "Jeb3vLu$^L8to=Dr{9YOxlIl'nw,XD")
]

COMMANDS = [
    "echo 'Тестовое подключение выполнено успешно'",
    "uname -a",
    "cat /etc/os-release",
    "df -h | head -n 5"
]


def main():
    tunnel = SSHTunnel()
    try:
        # Последовательное подключение через все хосты
        for host, port, user, pwd in JUMP_CHAIN:
            tunnel.connect(host, port, user, pwd)

        # Выполнение команд на конечном сервере
        results = tunnel.execute(COMMANDS)

        # Дополнительные действия с результатами
        with open('remote_results.txt', 'w') as f:
            for res in results:
                f.write(f"=== {res['command']} ===\n")
                f.write(res['output'] + "\n\n")

        print("\nРезультаты сохранены в remote_results.txt")

    except Exception as e:
        print(f"Ошибка: {str(e)}")
    finally:
        tunnel.close()


if __name__ == "__main__":
    main()