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