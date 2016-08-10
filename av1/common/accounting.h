/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
#ifndef AOM_ACCOUNTING_H_
#define AOM_ACCOUNTING_H_
#include <stdlib.h>

#define AOM_ACCOUNTING_HASH_SIZE (1021)

/* Max number of entries for symbol types in the dictionary (increase as
   necessary). */
#define MAX_SYMBOL_TYPES (256)

typedef struct {
  int16_t x;
  int16_t y;
  int8_t layer;
  int8_t level;
} AOMAccountingSymbolContext;

typedef struct {
  AOMAccountingSymbolContext context;
  uint32_t id;
  /** Number of bits in units of 1/8 bit. */
  int bits;
} AOMAccountingSymbol;

/** Dictionary for translating strings into id. */
typedef struct {
  char *(strs[MAX_SYMBOL_TYPES]);
  int num_strs;
} AOMAccountingDictionary;

typedef struct {
  /** All recorded symbols decoded. */
  AOMAccountingSymbol *syms;
  /** Number of symbols actually recorded. */
  int num_syms;
  /** Dictionary for translating strings into id. */
  AOMAccountingDictionary dictionary;
} AOMAccountingSymbols;

typedef struct {
  AOMAccountingSymbols syms;
  /** Size allocated for symbols (not all may be used). */
  int num_syms_allocated;
  int16_t hash_dictionary[AOM_ACCOUNTING_HASH_SIZE];
  AOMAccountingSymbolContext context;
  uint32_t last_tell;
} AOMAccounting;

void aom_accounting_init(AOMAccounting *accounting);
void aom_accounting_reset(AOMAccounting *accounting);
void aom_accounting_clear(AOMAccounting *accounting);
void aom_accounting_set_context(AOMAccounting *accounting, int16_t x, int16_t y);
int aom_accounting_dictionary_lookup(AOMAccounting *accounting, const char *str);
void aom_accounting_record(AOMAccounting *accounting, const char *str, int bits);
void aom_accounting_dump(AOMAccounting *accounting);
#endif  // AOM_ACCOUNTING_H_