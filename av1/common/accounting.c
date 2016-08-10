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

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "accounting.h"


static int aom_accounting_hash(const char *str) {
  uint32_t val;
  const unsigned char *ustr;
  val = 0;
  ustr = (const unsigned char*)str;
  /* This is about the worst hash one can design, but it should be good enough
     here. */
  while (*ustr) val += *ustr++;
  return val % AOM_ACCOUNTING_HASH_SIZE;
}

/* Dictionary lookup based on an open-addressing hash table. */
int aom_accounting_dictionary_lookup(AOMAccounting *accounting, const char *str) {
  int hash;
  AOMAccountingDictionary *dictionary;
  dictionary = &accounting->syms.dictionary;
  hash = aom_accounting_hash(str);
  while (accounting->hash_dictionary[hash] != -1) {
    if (strcmp(dictionary->strs[accounting->hash_dictionary[hash]], str) == 0) {
      return accounting->hash_dictionary[hash];
    }
    hash++;
    if (hash == AOM_ACCOUNTING_HASH_SIZE) hash = 0;
  }
  /* No match found. */
  assert(dictionary->num_strs + 1 < MAX_SYMBOL_TYPES);
  accounting->hash_dictionary[hash] = dictionary->num_strs;
  dictionary->strs[dictionary->num_strs] = malloc(strlen(str) + 1);
  strcpy(dictionary->strs[dictionary->num_strs], str);
  dictionary->num_strs++;
  return dictionary->num_strs - 1;
}

void aom_accounting_init(AOMAccounting *accounting) {
  int i;
  accounting->num_syms_allocated = 1000;
  accounting->syms.syms =
    malloc(sizeof(AOMAccountingSymbol) * accounting->num_syms_allocated);
  accounting->syms.dictionary.num_strs = 0;
  assert(AOM_ACCOUNTING_HASH_SIZE > 2 * MAX_SYMBOL_TYPES);
  for (i = 0; i < AOM_ACCOUNTING_HASH_SIZE; i++) accounting->hash_dictionary[i] = -1;
  aom_accounting_reset(accounting);
}

void aom_accounting_reset(AOMAccounting *accounting) {
  accounting->syms.num_syms = 0;
  accounting->context.x = -1;
  accounting->context.y = -1;
  accounting->context.level = -1;
  accounting->context.layer = -1;
  accounting->last_tell = 0;
}

void aom_accounting_clear(AOMAccounting *accounting) {
  int i;
  AOMAccountingDictionary *dictionary;
  free(accounting->syms.syms);
  dictionary = &accounting->syms.dictionary;
  for (i = 0; i < dictionary->num_strs; i++) {
    free(dictionary->strs[i]);
  }
}

void aom_accounting_set_context(AOMAccounting *accounting, int16_t x, int16_t y) {
  accounting->context.x = x;
  accounting->context.y = y;
}

void aom_accounting_record(AOMAccounting *accounting, const char *str, int bits) {
  AOMAccountingSymbol sym;
  sym.context = accounting->context;
  sym.bits = bits;
  sym.id = aom_accounting_dictionary_lookup(accounting, str);
  // assert(bits <= 255);
  assert(sym.id <= 255);
  if (accounting->syms.num_syms == accounting->num_syms_allocated) {
    accounting->num_syms_allocated *= 2;
    accounting->syms.syms = realloc(accounting->syms.syms,
      sizeof(AOMAccountingSymbol) * accounting->num_syms_allocated);
    assert(accounting->syms.syms != NULL);
  }
  accounting->syms.syms[accounting->syms.num_syms++] = sym;
}

void aom_accounting_dump(AOMAccounting *accounting) {
  int i;
  AOMAccountingSymbol *sym;
  printf("----- %d -----\n", accounting->syms.num_syms);
  for (i = 0; i < accounting->syms.num_syms; i++) {
    sym = &accounting->syms.syms[i];
    printf("%s %d, %d - %d\n", accounting->syms.dictionary.strs[sym->id], sym->context.x, sym->context.y, sym->bits);
  }
}