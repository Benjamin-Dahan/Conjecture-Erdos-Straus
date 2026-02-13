def est_premier(n):
    if n <= 1:
        return False
    if n == 2 or n == 3:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True

def check_conditions(n):
    if not est_premier(n):
        return False
    if n % 3 != 1 or n % 4 != 1 or (n % 7 != 4 and n % 7 != 2 and n % 7 != 1):
        return False
    if ((n + 3) / 4) % 2 == 0 or ((n + 3) / 4) % 3 == 0 or ((n + 3) / 4) % 5 == 0:
        return False
    return True


def main():
    for n in range(3362, 10000000000):
        if check_conditions(n):
            for c in range(3, 200001):
                if c % 4 == 3:
                    r = n % c
                    b = None

                    if r % 4 == 1:
                        b = 0
                    elif r % 4 == 2:
                        b = 1
                    elif r % 4 == 3:
                        b = 2
                    elif r % 4 == 0:
                        continue

                    sigma = (n + 4 * c - b * c - r) / (4 * c)

                    

                    csigma_minus_k = int(c * sigma - (n + c) / 4)

                    # Passer au c suivant si la condition sur congrq ou congrK n'est pas respectée
                    if (sigma * (n + c) / 4) % csigma_minus_k == 0:
                        print(f"n={n}, (n+c)/4={(n+c)/4}, c={c}, sigma={sigma}, ")
                        break
            else:
                print(f"Le programme s'arrête : aucun c n'a fonctionné pour n={n}")
                return

if __name__ == "__main__":
    main()
