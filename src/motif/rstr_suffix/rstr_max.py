#!/usr/bin/env python
# -*- coding: utf-8 -*-

from tools_karkkainen_sanders import direct_kark_sort
from array import array
import numpy as np

class Rstr_max :
  def __init__(self) :
    self.sequence = []
    self.array_str = []

  def add_str(self, seq) :
    self.sequence = seq
    self.array_str.append(seq)

  def step1and2(self):
    n = len(self.sequence)
    self.idxPos = np.arange(n)
    self.endAt = np.full(n, n)
    self.suffix_array = direct_kark_sort(self.sequence)
    rank = np.full(n, 0)
    self.lcp = np.full(n, 0)
    
    # figure out what this is doing
    for i in xrange(1,n):
      v = self.res[i]
      rank[v] = i
      print i,v

    l = 0
    for j in xrange(n):
      if(l > 0) :
        l -= 1
      i = rank[j]
      j2 = self.suffix_array[i-1]
      if i:
        while l + j < self.endAt[j] and l + j2 < self.endAt[j2] and self.sequence[j+l] == self.sequence[j2+l]:
          l += 1
        self.lcp[i-1] = l 
      else:
        l = 0

  def step3_rstr(self) :
    prev_len = idx = 0
    results = {}
    len_lcp = len(self.lcp) -1

    self.top = 0
    self.lst_max = []

    if len(self.suffix_array) == 0 :
      return {}

    pos1 = self.suffix_array[0]
    for idx in xrange(len_lcp):
      current_len = self.lcp[idx]
      pos2 = self.res[idx+1]
      end = max(pos1, pos2) + current_len
      n = prev_len - current_len
      if n < 0 :
        self.lst_max.append([-n, idx, end])
        self.top += -n
      elif n > 0:
        self.removeMany(results, n, idx)
      elif self.top > 0 and end > self.lst_max[-1][-1] :
        self.lst_max[-1][-1] = end
      prev_len = current_len
      pos1 = pos2
      
    if(self.top > 0) :
      self.removeMany(results, self.top, idx+1)

    return results

  def removeMany(self, results, m, idxEnd):
      prevStart = -1
      while m > 0:
        n, idxStart, maxEnd = self.lst_max.pop()
        if prevStart != idxStart:
          id_ = (maxEnd, idxEnd-idxStart+1)
          if id_ not in results or results[id_][0] < self.top:
              results[id_] = (self.top, idxStart)
          prevStart = idxStart
        m -= n
        self.top -= n
      if m < 0:
        self.lst_max.append([-m, idxStart, maxEnd-n-m])
        self.top -= m    
    
  def step1_sort_suffix(self) :
    char_frontier = chr(2)
    self.global_suffix = char_frontier.join(self.array_str)
    nbChars = len(self.global_suffix)
    init = [-1]*nbChars
    self.idxString = array('i', init)
    self.idxPos = array('i', init)
    self.endAt = array('i', init)
    # self.idxPos = np.arange(nbChars)
    # self.endAt = np.full(nbChars, nbChars)
    k = idx = 0
    for mot in self.array_str :
      last = k + len(mot)
      for p in xrange(len(mot)) :
        self.idxString[k] = idx
        self.idxPos[k] = p
        self.endAt[k] = last
        k += 1
      idx += 1
      k += 1

    self.res = direct_kark_sort(self.global_suffix)

  def step2_lcp(self) :
    n = len(self.res)
    init = [0]*n
    rank = array('i', init)
    LCP = array('i', init)
    print rank, LCP
    
    s = self.global_suffix
    suffix_array = self.res
    endAt = self.endAt

    for i in xrange(len(self.array_str),n):
      v = self.res[i]
      rank[v] = i
      print i,v
    print self.res
    print rank
    l = 0
    for j in xrange(n):
      if(l > 0) :
        l -= 1
      i = rank[j]
      j2 = suffix_array[i-1]
      if i:
        while l + j < endAt[j] and l + j2 < endAt[j2] and s[j+l] == s[j2+l]:
          l += 1
        LCP[i-1] = l 
      else:
        l = 0
    self.lcp = LCP

  def step3_rstr(self) :
    prev_len = 0
    idx = 0
    results = {}
    len_lcp = len(self.lcp) -1

    class Stack:
      pass
    stack = Stack()
    stack._top = 0
    stack.lst_max = []

    if len(self.res) == 0 :
      return {}

    pos1 = self.res[0]
    for idx in xrange(len_lcp):
      current_len = self.lcp[idx]
      pos2 = self.res[idx+1]
      end_ = max(pos1, pos2) + current_len
      n = prev_len - current_len
      if n < 0 :
        stack.lst_max.append([-n, idx, end_])
        stack._top += -n
      elif n > 0:
        self.removeMany(stack, results, n, idx)
      elif stack._top > 0 and end_ > stack.lst_max[-1][-1] :
        stack.lst_max[-1][-1] = end_

      prev_len = current_len
      pos1 = pos2
      
    if(stack._top > 0) :
      self.removeMany(stack, results, stack._top, idx+1)

    return results

  def removeMany(self, stack, results, m, idxEnd):
      prevStart = -1
      while m > 0:
        n, idxStart, maxEnd = stack.lst_max.pop()
        if prevStart != idxStart:
          id_ = (maxEnd, idxEnd-idxStart+1)
          if id_ not in results or results[id_][0] < stack._top:
              results[id_] = (stack._top,idxStart)
          prevStart = idxStart
        m -= n
        stack._top -= n
      if m < 0:
        stack.lst_max.append([-m, idxStart, maxEnd-n-m])
        stack._top -= m    
    
  def go(self) :
    self.step1_sort_suffix()
    self.step2_lcp()
    r = self.step3_rstr()
    return r
  

def GetMotifs(sequence):
  '''
  Given a sequence, return a list of:
  # (length), [<start_indices>]
  '''
  unic = unicode(sequence, 'utf-8', 'replace')
  rstr = Rstr_max()
  rstr.add_str(unic)
  r = rstr.go()
  result = [] # (word_index_start, word_index_end) 
  for (_, nb), (l, start) in r.iteritems():
    occurrences = [rstr.idxPos[rstr.res[o]] for o in range(start, start+nb)]
    result.append((l, occurrences))
  return result


if (__name__ == '__main__') :
  str1 = 'ababtototiti'
  print GetMotifs(str1)


