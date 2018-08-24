"""
@summary: This module contains classes and functions for testing the provided
             data readers
@note: Uses pytest style testing
"""
import pytest

# .............................................................................
class Test_create_sequence_list_from_dict(object):
   def test_invalid_dict(self):
      pass
   
   def test_valid_dict(self):
      pass

# .............................................................................
class Test_read_csv_alignment_flo(object):
   def test_file_invalid(self):
      pass
   
   def test_file_valid(self, data_files):
      print data_files.VALID_CSV_ALIGNMENT
      assert False
      
      pass
   
   def test_stringio_invalid(self):
      pass
   
   def test_stringio_valid(self):
      pass

# .............................................................................
class Test_read_json_alignment_flo(object):
   def test_file_invalid(self):
      pass
   
   def test_file_valid(self):
      pass
   
   def test_stringio_invalid(self):
      pass
   
   def test_stringio_valid(self):
      pass

# .............................................................................
class Test_read_phylip_alignment_flo(object):
   def test_file_invalid(self):
      pass
   
   def test_file_valid(self):
      pass
   
   def test_stringio_invalid(self):
      pass
   
   def test_stringio_valid(self):
      pass

# .............................................................................
class Test_read_phylip_continuous_alignment_flo(object):
   def test_file_invalid(self):
      pass
   
   def test_file_valid(self):
      pass
   
   def test_stringio_invalid(self):
      pass
   
   def test_stringio_valid(self):
      pass

# .............................................................................
class Test_read_table_continuous_alignment_flo(object):
   def test_file_invalid(self):
      pass
   
   def test_file_valid(self):
      pass
   
   def test_stringio_invalid(self):
      pass
   
   def test_stringio_valid(self):
      pass

# .............................................................................

# content of test_class.py
class TestClass(object):
    def test_one(self):
        x = "this"
        assert 'h' in x

    def test_two(self):
        x = "hello"
        assert hasattr(x, 'check')

def test_cj(data_files):
   print data_files.VALID_CSV_ALIGNMENT
   assert True
   assert False
   
   