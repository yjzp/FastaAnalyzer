"""
Расширенные тесты для класса FastaReader с использованием pytest.
"""

import pytest
import sys
import os
import tempfile

# Добавляем путь к модулю
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from fasta_parser.fasta_reader import FastaReader
from fasta_parser.seq import Seq
from fasta_parser.exceptions import FastaFormatError, FileNotFoundError


class TestFastaReaderExtended:
    """Расширенные тесты для класса FastaReader."""

    @pytest.fixture
    def test_data_dir(self):
        """Фикстура для получения пути к тестовым данным."""
        return os.path.join(os.path.dirname(__file__), 'test_data')

    @pytest.fixture
    def sample_files(self, test_data_dir):
        """Фикстура с путями к тестовым файлам."""
        return {
            'dna': os.path.join(test_data_dir, 'sample_dna.fasta'),
            'protein': os.path.join(test_data_dir, 'sample_protein.fasta'),
            'rna': os.path.join(test_data_dir, 'sample_rna.fasta'),
            'large': os.path.join(test_data_dir, 'large_test.fasta'),
            'invalid': os.path.join(test_data_dir, 'invalid.txt'),
            'error': os.path.join(test_data_dir, 'error_test.fasta')
        }

    def test_file_validation(self, sample_files):
        """Тест валидации файлов."""
        # Корректные FASTA файлы
        for file_type in ['dna', 'protein', 'rna']:
            if os.path.exists(sample_files[file_type]):
                reader = FastaReader(sample_files[file_type])
                assert reader.is_fasta_format()

        # Некорректный файл
        if os.path.exists(sample_files['invalid']):
            reader = FastaReader(sample_files['invalid'])
            assert not reader.is_fasta_format()

    def test_nonexistent_file(self):
        """Тест с несуществующим файлом."""
        with pytest.raises(FileNotFoundError):
            FastaReader("nonexistent_file.fasta")

    def test_sequence_counting(self, sample_files):
        """Тест подсчета последовательностей."""
        if os.path.exists(sample_files['dna']):
            reader = FastaReader(sample_files['dna'])
            count = reader.count_sequences()
            assert count > 0
            assert count == len(reader)

    def test_sequence_reading(self, sample_files):
        """Тест чтения последовательностей."""
        if os.path.exists(sample_files['dna']):
            reader = FastaReader(sample_files['dna'])
            sequences = list(reader.read_sequences())

            assert len(sequences) > 0
            for seq in sequences:
                assert isinstance(seq, Seq)
                assert len(seq) > 0
                assert seq.header.strip() != ""

    def test_iteration(self, sample_files):
        """Тест итерации по файлу."""
        if os.path.exists(sample_files['dna']):
            reader = FastaReader(sample_files['dna'])

            count = 0
            for seq in reader:
                assert isinstance(seq, Seq)
                count += 1

            assert count == len(reader)

    def test_sequence_search(self, sample_files):
        """Тест поиска последовательностей."""
        if os.path.exists(sample_files['dna']):
            reader = FastaReader(sample_files['dna'])

            # Получаем первую последовательность для поиска
            first_seq = next(reader.read_sequences())
            search_term = first_seq.header.split()[0]  # Первое слово заголовка

            found_seq = reader.get_sequence_by_id(search_term)
            assert found_seq is not None
            assert search_term in found_seq.header

    def test_type_filtering(self, sample_files):
        """Тест фильтрации по типу."""
        if os.path.exists(sample_files['dna']):
            reader = FastaReader(sample_files['dna'])
            dna_sequences = list(reader.get_sequences_by_type("DNA"))

            assert len(dna_sequences) > 0
            for seq in dna_sequences:
                assert seq.get_alphabet_type() == "DNA"

    def test_file_statistics(self, sample_files):
        """Тест получения статистики файла."""
        if os.path.exists(sample_files['dna']):
            reader = FastaReader(sample_files['dna'])
            stats = reader.get_file_stats()

            assert 'total_sequences' in stats
            assert 'total_length' in stats
            assert 'sequence_types' in stats
            assert 'length_stats' in stats

            assert stats['total_sequences'] > 0
            assert stats['total_length'] > 0
            assert stats['length_stats']['avg_length'] > 0

    def test_large_file_processing(self, sample_files):
        """Тест обработки большого файла."""
        if os.path.exists(sample_files['large']):
            reader = FastaReader(sample_files['large'])

            # Должен работать как генератор
            count = 0
            total_length = 0

            for seq in reader:
                count += 1
                total_length += len(seq)
                # Проверяем только первые 10 для экономии времени
                if count >= 10:
                    break

            assert count > 0
            assert total_length > 0

    def test_filtered_output(self, sample_files):
        """Тест записи отфильтрованного вывода."""
        if os.path.exists(sample_files['dna']):
            reader = FastaReader(sample_files['dna'])

            # Создаем временный файл
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_file:
                temp_filename = temp_file.name

            try:
                # Фильтр: только длинные последовательности
                def long_sequences_filter(seq, min_length=50):
                    return len(seq) >= min_length

                written_count = reader.write_filtered_fasta(
                    temp_filename,
                    filter_func=long_sequences_filter,
                    min_length=50
                )

                assert written_count >= 0

                # Проверяем результат
                if written_count > 0:
                    filtered_reader = FastaReader(temp_filename)
                    assert filtered_reader.is_fasta_format()

                    filtered_sequences = list(filtered_reader.read_sequences())
                    assert len(filtered_sequences) == written_count

                    for seq in filtered_sequences:
                        assert len(seq) >= 50

            finally:
                # Удаляем временный файл
                if os.path.exists(temp_filename):
                    os.unlink(temp_filename)

    def test_error_handling(self, sample_files):
        """Тест обработки ошибок."""
        # Тест с невалидным FASTA файлом
        if os.path.exists(sample_files['invalid']):
            reader = FastaReader(sample_files['invalid'])

            with pytest.raises(FastaFormatError):
                list(reader.read_sequences())

    def test_empty_sequences(self):
        """Тест с пустыми последовательностями."""
        empty_fasta_content = """>empty_seq1

>empty_seq2

>normal_seq
ATGC
"""

        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_file:
            temp_file.write(empty_fasta_content)
            temp_filename = temp_file.name

        try:
            reader = FastaReader(temp_filename)
            sequences = list(reader.read_sequences(validate=False))  # Отключаем валидацию

            # Должно быть 3 последовательности
            assert len(sequences) == 3

            # Первые две пустые
            assert len(sequences[0]) == 0
            assert len(sequences[1]) == 0

            # Третья нормальная
            assert len(sequences[2]) == 4
            assert sequences[2].sequence == "ATGC"

        finally:
            if os.path.exists(temp_filename):
                os.unlink(temp_filename)

    def test_different_encodings(self, sample_files):
        """Тест с различными кодировками."""
        if os.path.exists(sample_files['dna']):
            # Тест с UTF-8 (по умолчанию)
            reader_utf8 = FastaReader(sample_files['dna'], encoding='utf-8')
            sequences_utf8 = list(reader_utf8.read_sequences())

            assert len(sequences_utf8) > 0

    def test_memory_efficiency(self, sample_files):
        """Тест эффективности использования памяти."""
        if os.path.exists(sample_files['large']):
            reader = FastaReader(sample_files['large'])

            # Генератор не должен загружать все в память сразу
            seq_iterator = reader.read_sequences()

            # Получаем только первую последовательность
            first_seq = next(seq_iterator)
            assert isinstance(first_seq, Seq)

            # Можем продолжить итерацию без проблем
            try:
                second_seq = next(seq_iterator)
                assert isinstance(second_seq, Seq)
            except StopIteration:
                # Файл может содержать только одну последовательность
                pass

    @pytest.mark.parametrize("seq_type,expected_min_count", [
        ("DNA", 1),
        ("RNA", 0),
        ("PROTEIN", 0),
    ])
    def test_sequence_type_counting(self, sample_files, seq_type, expected_min_count):
        """Параметризованный тест подсчета типов последовательностей."""
        if os.path.exists(sample_files['dna']):
            reader = FastaReader(sample_files['dna'])
            sequences = list(reader.get_sequences_by_type(seq_type))
            assert len(sequences) >= expected_min_count


if __name__ == '__main__':
    pytest.main([__file__])
