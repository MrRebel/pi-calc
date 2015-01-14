with open("pi.txt") as f1, open("realpi.txt") as f2:
    pi = f1.read()
    realpi = f2.read()

    realpi = ''.join((i for i in realpi if i in '0123456789'))
   
i = 0
try:
    while pi[i] == realpi[i]:
        i += 1
except IndexError:
    pass
print(i, "digits out of", min((len(pi), len(realpi))) ,"verified correct!")
