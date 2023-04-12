#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <cstdlib>

using namespace std;

ifstream in ("date.in");
ofstream out("Evolutie.txt");

int n, a, b, c[3], precizie, nr_etape, len, ind, power;
long double maxx, p_recombinare, p_mutatie, minn;
vector <long long> chromosomes, new_gen;
vector <long double> vals, probs, intervale;
vector <int> viz;
vector <int> recombinare, aux;

// functie de afisare a numerelor in baza 2
void bin(long long nr) {
    long long i;
    for (i = 1 << (len - 1); i > 0; i = i / 2) {
        if((nr & i) != 0) {
            out << "1";
        }
        else {
            out << "0";
        }
    }
}

// functie de rotunjire cu pow zecimale
long double round_nr (long double nr, int pow) {
    int value = (int)(nr * (long double)pow);
    return (long double)value / (long double)pow;
}

long double functie (long double x) {
    return ((long double)c[0] * x * x + c[1] * x + c[2]);
}

int cautare_binara (long double nr) {
    int mijl, st = 0, dr = n;
    while (st <= dr) {
        mijl = (st + dr);
        if (intervale[mijl] <= nr) {
            st = mijl + 1;
        } else if (intervale[mijl] > nr) {
            dr = mijl - 1;
        }
    }
    return dr;
}

int main() {

    cout << "Introduceti dimensiunea populatiei: "; in >> n;
    cout << "Introduceti domeniul de definitie al functiei: "; in >> a >> b;
    cout << "Introduceti coeficientii functiei: "; in >> c[0] >> c[1] >> c[2];
    cout << "Introduceti precizia cu care se lucreaza: "; in >> precizie;
    cout << "Introduceti probabilitatea de recombinare: "; in >> p_recombinare;
    cout << "Introduceti probabilitatea de mutatie: "; in >> p_mutatie;
    cout << "Introduceti numarul de etape al algoritmului: "; in >> nr_etape;
    cout<<'\n';

    power = pow(10, precizie);
    len = ceil(log2((b - a) * power)); // lungimea cromozomului
    long double cst = ((long double)b - a) / (pow(2, len) - 1);


    // la primul pas alegem indivizii populatiei random
    for (int i = 0; i < n; i++) {
        chromosomes.push_back(rand() % (1 << len));
    }


    for (int gen = 0; gen < nr_etape; gen++) {
        if (gen == 0) {
            out << "Populatia initiala:\n";
        }
        long double total_fitness = 0;
        vals.clear();
        probs.clear();
        intervale.clear();
        viz.clear();
        recombinare.clear();
        //maxx = 0;
        long long ch_fit;
        for (int i = 0; i < n; i++) {
            long long y = chromosomes[i];
            long double val_y = cst * y + a;  // valoarea codificata din interval - translatie liniara
            vals.push_back(val_y);
            //long double fct = functie(round_nr(vals[i], power));

            if (gen == 0) {
                out << i + 1 << ": ";
                bin(y);
                out << " x= " << round_nr(vals[i], power) << " f= " << functie(round_nr(vals[i], power)) << '\n';
            }
            total_fitness += functie(vals[i]);
            if (i == 0) {
                maxx = functie(round_nr(vals[i], power));
                ind = i;
            }
            else if (functie(round_nr(vals[i], power)) > maxx) {
                maxx = functie(round_nr(vals[i], power));
                ind = i;
            }
        }

        long double max_fit = maxx; // cel mai mare scor fitness
        ch_fit = chromosomes[ind]; // cel mai fit individ din generatie

        if (gen == 0) {
            out << "\nProbabilitati selectie:\n";
        }
        for (int i = 0; i < n; i++) {
            long double prob_i = functie(vals[i]) / total_fitness;
            probs.push_back(prob_i);
            if (gen == 0) {
                out << "cromozom " << i + 1 << " probabilitate " << prob_i << '\n';
            }
        }

        if (gen == 0) {
            out << "\nIntervale probabilitati selectie:\n";
        }
        intervale.push_back(0);
        long double sum_interv = probs[0];
        intervale.push_back(sum_interv);

        if (gen == 0) {
            out << 0 << " " << sum_interv << " ";
        }
        for (int i = 1; i < n; i++) {
            sum_interv += probs[i];
            intervale.push_back(sum_interv);
            if (gen == 0) {
                out << sum_interv << ' ';
            }
        }

        // selectia proportionala - metoda ruletei
        if (gen == 0) out << '\n';
        for (int i = 0; i < n; i++) {
            long double u = static_cast <long double> (rand()) / static_cast <long double> (RAND_MAX);
            int cromozom = cautare_binara(u);
            if (gen == 0) {
                out << "u= " << u << " selectam cromozomul " << cromozom + 1 << '\n';
            }
            viz.push_back(cromozom);
        }


        if (gen == 0) {
            out << "\nDupa selectie:\n";
        }
        for (int i = 0; i < n; i++) {
            if (gen == 0) {
                out << i + 1 << ": ";
                bin(chromosomes[viz[i]]);
                out << " x= " << vals[viz[i]] << " f= " << functie(vals[viz[i]]) << '\n';
            }
            new_gen.push_back(chromosomes[viz[i]]);
        }
        chromosomes = new_gen;

        // selectam cromozomii ce vor fi recombinati
        if (gen == 0) {
            out << "\nProbabilitatea de incrucisare: " << p_recombinare << '\n';
        }
        for (int i = 0; i < n; i++) {
            long double u = static_cast <long double> (rand()) / static_cast <long double> (RAND_MAX);
            if (u < p_recombinare) {
                recombinare.push_back(i);
                if (gen == 0) {
                    out << i + 1 << ": ";
                    bin(chromosomes[i]);
                    out << " u= " << u << " < " << p_recombinare << " participa\n";
                }
            }
            else {
                if (gen == 0) {
                    out << i + 1 << ": ";
                    bin(chromosomes[i]);
                    out << " u= " << u << "\n";
                }
            }
        }

        // recombinam cate 2 cromozomi random
        while (recombinare.size() > 1) {
            int x = rand() % recombinare.size();
            int y = rand() % recombinare.size();
            if (x == y) continue;
            if (gen == 0) {
                out << "Recombinare dintre cromozomul " << recombinare[x] + 1 << " cu cromozomul " << recombinare[y] + 1 << '\n';
            }
            int linie = rand() % len;
            if (gen == 0) {
                bin(chromosomes[recombinare[x]]);
                out << ' ';
                bin(chromosomes[recombinare[y]]);
                out << " punct " << linie << "\nRezultat ";

            }
            long long masca = (1 << (len - linie + 1)) -1;
            long long masca_x = chromosomes[recombinare[x]] & masca;
            long long masca_y = chromosomes[recombinare[y]] & masca;
            chromosomes[recombinare[x]] >>= (len - linie);
            chromosomes[recombinare[x]] <<= (len - linie);
            chromosomes[recombinare[x]] |= masca_y;
            chromosomes[recombinare[y]] >>= (len - linie);
            chromosomes[recombinare[y]] <<= (len - linie);
            chromosomes[recombinare[y]] |= masca_x;

            if (gen == 0) {
                bin(chromosomes[recombinare[x]]);
                out << ' ';
                bin(chromosomes[recombinare[y]]);
                out << '\n';
            }

            // elimin cromozomii abia recombinati din vectorul de recombinare
            aux.clear();
            for (int i = 0; i < recombinare.size(); i++) {
                if (i != x && i != y) {
                    aux.push_back(recombinare[i]);
                }
            }
            recombinare = aux;
        }

        if (gen == 0) {
            out << "\nDupa recombinare:\n";
        }

        for (int i = 0; i < n; i++) {

            // recalculez valoarea codificata din interval
            long long y = chromosomes[i];
            long double val_y = cst * y + a;
            vals[i] = val_y;

            if (gen == 0) {
                out << i + 1 << ": ";
                bin(y);
                out << " x= " << vals[i] << " f= " << round_nr(functie(vals[i]), power) << '\n';
            }
        }

        if (gen == 0) {
            out << "\nProbabilitatea de mutatie " << p_mutatie << '\n';
            out << "Au fost modificati cromozomii:\n";
        }
        for (int i = 0; i < n; i++) {
            long double u = static_cast <long double> (rand()) / static_cast <long double> (RAND_MAX);
            if (u < p_mutatie) {
                int poz = rand() % len;
                long long masca = (1 << (len - poz + 1)) - 1;
                masca &= chromosomes[i];
                chromosomes[i] >>= (len - poz);
                if (chromosomes[i] % 2 == 0) chromosomes[i]++;
                else chromosomes[i]--;
                chromosomes[i] <<= (len - poz);
                chromosomes[i] |= masca;

                if (gen == 0) {
                    out << i + 1 << '\n';
                }
            }
        }

        if (gen == 0) {
            out << "\nDupa mutatie:\n";
        }

        minn = 0;
        maxx = 0;
        for (int i = 0; i < n; i++) {

            // recalculez valoarea codificata din interval
            long long y = chromosomes[i];
            long double val_y = cst * y + a;
            vals[i] = val_y;
            long double fct = functie(vals[i]);
            if (gen == 0) {
                out << i + 1 << ": ";
                bin(y);
                out << " x= " << vals[i] << " f= " << round_nr(fct, power) << '\n';
            }
            if (i == 0) {
                maxx = fct;
                minn = fct;
            }
            else {
                // retin maximul de fitness
                if (fct > maxx) {
                    maxx = fct;
                }
                // retin cel mai nefit cromozom
                if (fct < minn) {
                    minn = fct;
                    ind = i;
                }
            }
        }

        // selectia elitista
        // schimb cel mai nefit cromozom cu cel mai fit din prima generatie inainte de schimbari
        chromosomes[ind] = ch_fit;


        if (max_fit > maxx) {
            maxx = max_fit;
        }

        if (gen == 0) {
            out << "\nEvolutia maximului:\n";
        }

        out << maxx << '\n';

    }

    return 0;
}
