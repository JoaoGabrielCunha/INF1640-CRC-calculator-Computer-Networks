
/*
 * CRC + LFSR demonstrator (C version)
 * Polinômio: g(x) = x^6 + x^4 + x^3 + x + 1  (0b1011011)
 *
 * O programa:
 *  (1) Faz a divisão em GF(2) mostrando os passos (quociente e resto).
 *  (2) Calcula a mensagem transmitida (codeword) e verifica o resto = 0.
 *  (3) Gera a tabela de evolução do LFSR (32 bits + 6 zeros) e compara FCS.
 * Também grava toda a saída em um arquivo texto além do stdout.
 *
 * Compile:  gcc -std=c11 -O2 -Wall -Wextra -o crc_lfsr crc_lfsr.c
 */

#include <stdio.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

/* ===================== util: logger duplo (stdout + arquivo) ===================== */
typedef struct {
    FILE *fp;
} Logger;

static void lprint(Logger *L, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    vprintf(fmt, ap);
    va_end(ap);

    if (L && L->fp) {
        va_start(ap, fmt);
        vfprintf(L->fp, fmt, ap);
        va_end(ap);
    }
}

static void print_repeat(Logger *L, char c, int n) {
    for (int i = 0; i < n; ++i) lprint(L, "%c", c);
}

static int bitlen_u64(uint64_t x) {
    if (x == 0) return 0;
    int n = 0;
    while (x) { n++; x >>= 1; }
    return n;
}

/* Constrói string "0b" + <width bits de x> (com zeros à esquerda). Retorna malloc'd. */
static char* bits_str(uint64_t x, int width) {
    int len = 2 + (width > 0 ? width : 1);
    char *s = (char*)malloc((len + 1) * sizeof(char));
    if (!s) return NULL;
    s[0] = '0'; s[1] = 'b';
    if (width <= 0) {
        s[2] = '0'; s[3] = '\0';
        return s;
    }
    for (int i = 0; i < width; ++i) {
        int bit = (int)((x >> (width - 1 - i)) & 1ULL);
        s[2 + i] = (char)('0' + bit);
    }
    s[len] = '\0';
    return s;
}

/* ===================== (1) Divisão em GF(2) com passos ===================== */
static void divide_mod2_show(uint64_t dividendo, uint64_t divisor,
                             uint64_t *q_out, uint64_t *r_out,
                             Logger *L, int verbose)
{
    if (divisor == 0) {
        lprint(L, "Erro: divisor não pode ser zero.\n");
        exit(1);
    }
    uint64_t quoc = 0;
    int k = bitlen_u64(dividendo);           /* bits do dividendo */
    int r = bitlen_u64(divisor);             /* bits do divisor   */
    int m = (r > 0 ? r - 1 : 0);             /* grau do polinômio */
    int aux = k - r;                         

    int binw = 2 + r;                        /* largura '0b' + r bits */

    if (aux < 0) {
        uint64_t resto = dividendo;
        if (verbose) {
            lprint(L, "Divisão módulo 2 (dividendo menor que divisor)\n");
            char *sd = bits_str(dividendo, k);
            char *sr = bits_str(divisor,   r);
            lprint(L, "%s |__ %s\n", sd, sr);
            free(sd); free(sr);
            lprint(L, "Quociente: 0b0\n");
            char *sre = bits_str(resto, r);
            lprint(L, "Resto:     %s\n", sre);
            free(sre);
        }
        if (q_out) *q_out = 0;
        if (r_out) *r_out = resto;
        return;
    }

    uint64_t resto = dividendo >> aux;

    if (verbose) {
        lprint(L, "Divisão módulo 2\n");
        char *sd = bits_str(dividendo, k);
        char *sr = bits_str(divisor,   r);
        lprint(L, "%s |__ %s\n", sd, sr);
        free(sd); free(sr);
    }

    while (aux > -1) {
        int nbits_resto = bitlen_u64(resto);
        quoc <<= 1;
        if (nbits_resto == r) {
            quoc |= 1;
            resto ^= divisor;
            if (verbose) {
                print_repeat(L, ' ', k - r - aux);
                char *s = bits_str(divisor, r);
                lprint(L, "%*s", binw, s);
                free(s);
                print_repeat(L, '|', aux);
                lprint(L, "\n");
            }
        } else {
            if (verbose) {
                print_repeat(L, ' ', k - r - aux);
                char *zero_line = bits_str(0, r);
                lprint(L, "%*s", binw, zero_line);
                free(zero_line);
                print_repeat(L, '|', aux);
                lprint(L, "\n");
            }
        }

        if (verbose) {
            print_repeat(L, ' ', k - r - aux);
            print_repeat(L, '-', binw);
            print_repeat(L, '|', aux);
            lprint(L, "\n");
        }

        aux -= 1;
        if (aux > -1) {
            uint64_t next_bit = (dividendo >> aux) & 1ULL;
            resto = (resto << 1) | next_bit;
        }

        if (verbose) {
            print_repeat(L, ' ', k - r - aux);
            char *sr = bits_str(resto, r);
            lprint(L, "%*s", binw, sr);
            free(sr);
            print_repeat(L, '|', aux);
            lprint(L, "\n");
        }
    }

    if (verbose) {
        lprint(L, "\nQuociente: ");
        char *sq = bits_str(quoc, r);
        lprint(L, "%s\n", sq); free(sq);
        lprint(L, "Resto: ");
        char *sr = bits_str(resto, r);
        lprint(L, "%s\n", sr); free(sr);
    }

    if (q_out) *q_out = quoc;
    if (r_out) *r_out = resto;
}

/* ===================== (1b) Constrói codeword e FCS ===================== */
static void make_crc_transmission(uint64_t mensagem, uint64_t polinomio,
                                  uint64_t *codeword_out, uint64_t *fcs_out,
                                  Logger *L, int verbose)
{
    int m = bitlen_u64(polinomio) - 1;
    uint64_t shifted = mensagem << m;
    uint64_t quo=0, rem=0;
    divide_mod2_show(shifted, polinomio, &quo, &rem, L, verbose);
    uint64_t codeword = shifted ^ rem;
    if (codeword_out) *codeword_out = codeword;
    if (fcs_out) *fcs_out = rem;
}

/* ===================== (2/3) LFSR: shift-in + XOR (MSB-old) ===================== */

static uint64_t trace_lfsr_crc(uint64_t mensagem, int msg_width,
                               uint64_t polinomio, Logger *L, int verbose)
{
    int r = bitlen_u64(polinomio);
    int m = r - 1;
    uint64_t mask_m = (m > 0) ? ((1ULL << m) - 1ULL) : 0;
    uint64_t poly_lo = polinomio & mask_m;
    uint64_t reg = 0;

    /* helper: imprime registrador como m bits */
    char *(*reg_str)(uint64_t, int) = bits_str;

    if (verbose) {
        lprint(L, "passo | i | msb(old) |  r[m-1]..r[0]      ->   r'[m-1]..r'[0]\n");
    }

    for (int step = 0; step < msg_width; ++step) {
        int i = (int)((mensagem >> (msg_width - 1 - step)) & 1ULL);
        int msb_old = (m>0) ? (int)((reg >> (m-1)) & 1ULL) : 0;
        uint64_t before = reg;
        reg = ((reg << 1) | (uint64_t)i) & mask_m;
        if (msb_old) reg ^= poly_lo;

        if (verbose) {
            char *b1 = reg_str(before, m);
            char *b2 = reg_str(reg,    m);
            lprint(L, "%5d | %d |     %d     |  %-16s ->   %s\n", step, i, msb_old, b1+2, b2+2);
            free(b1); free(b2);
        }
    }

    for (int z = 0; z < m; ++z) {
        int msb_old = (m>0) ? (int)((reg >> (m-1)) & 1ULL) : 0;
        uint64_t before = reg;
        reg = (reg << 1) & mask_m;
        if (msb_old) reg ^= poly_lo;

        if (verbose) {
            char *b1 = reg_str(before, m);
            char *b2 = reg_str(reg,    m);
            lprint(L, "%5d | 0 |     %d     |  %-16s ->   %s\n", msg_width+z, msb_old, b1+2, b2+2);
            free(b1); free(b2);
        }
    }

    return reg; /* FCS */
}


static void print_bits(Logger *L, const char *label, uint64_t x, int width) {
    char *s = bits_str(x, width);
    lprint(L, "%s%s\n", label, s);
    free(s);
}


int main(void) {
    /* Dados do enunciado */
    uint64_t mensagem  = 0b10001000100010001000000110000001ULL; /* 32 bits */
    uint64_t polinomio = 0b1011011ULL;                          /* x^6 + x^4 + x^3 + x + 1 */

    Logger logger = {0};
    logger.fp = fopen("resultado_crc.txt", "w");
    if (!logger.fp) {
        fprintf(stderr, "Aviso: não consegui abrir resultado_crc.txt para escrita.\n");
    }

    int m = bitlen_u64(polinomio) - 1;
    int msgw = 32;

    lprint(&logger, "\n=== ITEM 1: CRC por divisão em módulo 2 (com passos) ===\n\n");
    uint64_t codeword = 0, fcs_div = 0;
    make_crc_transmission(mensagem, polinomio, &codeword, &fcs_div, &logger, 1);

    lprint(&logger, "\n");
    print_bits(&logger, "FCS (divisão): ", fcs_div, m);
    {
        int cw_w = msgw + m;
        char *cw = bits_str(codeword, cw_w);
        lprint(&logger, "Mensagem transmitida (codeword): %s\n\n", cw);
        free(cw);
    }

    lprint(&logger, "Verificação na recepção (codeword ÷ polinômio):\n");
    {
        uint64_t q=0, r=0;
        divide_mod2_show(codeword, polinomio, &q, &r, &logger, 1);
        lprint(&logger, "\n%s\n", (r == 0) ? "Transmissão com sucesso!" : "Falha na transmissão.");
    }

    lprint(&logger, "\n=== ITEM 2 e 3: LFSR simplificado + tabela de evolução ===\n\n");
    uint64_t fcs_lfsr = trace_lfsr_crc(mensagem, msgw, polinomio, &logger, 1);

    lprint(&logger, "\n");
    print_bits(&logger, "FCS (LFSR):     ", fcs_lfsr, m);
    lprint(&logger, "Comparação:     %s\n\n", (fcs_lfsr == fcs_div) ? "OK" : "DIVERGE");

    if (logger.fp) fclose(logger.fp);
    return 0;
}
