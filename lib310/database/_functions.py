import random
import string

def number_to_abbrevation(num):
    label = ['', 'K', 'M', 'B', 'T']
    i = 0
    while num > 1000:
        num /= 1000
        i += 1
    return f'{num:.1f} {label[i]}'


def size_to_abbrevation(num):
    label = ['B', 'KB', 'MB', 'GB', 'TB']
    i = 0
    while num > 1024:
        num /= 1024
        i += 1
    return f'{num:.1f} {label[i]}'


def get_random_string(length):
    return ''.join(random.choice(string.ascii_uppercase) for i in range(length))