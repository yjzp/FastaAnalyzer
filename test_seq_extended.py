"""
Расширенные тесты для класса Seq с использованием pytest.
"""

import pytest
import sys
import os

# Добавляем путь к модулю
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from fasta_parser.seq import Seq
from fasta_parser.exceptions import SequenceError


class TestSeqExtended:
    """Расширенные тесты для класса Seq."""

    def test_dna_sequences(self):
        """Тест различных ДНК последовательностей."""
        # Обычная ДНК
        dna = Seq("ATGCGATCG", "Normal DNA")
        assert dna.get_alphabet_type() == "DNA"
        assert len(dna) == 9

        # ДНК с неоднозначными символами
        ambiguous_dna = Seq("ATGCRYSWKMBD", "Ambiguous DNA")
        assert ambiguous_dna.get_alphabet_type() == "DNA"

        # Длинная ДНК
        long_dna = Seq("A" * 1000 + "T" * 1000 + "G" * 1000 + "C" * 1000, "Long DNA")
        assert len(long_dna) == 4000
        assert long_dna.get_gc_content() == 50.0

    def test_rna_sequences(self):
        """Тест РНК последовательностей."""
        rna = Seq("AUGCGAUCG", "Normal RNA")
        assert rna.get_alphabet_type() == "RNA"

        # РНК не должна содержать T
        mixed_rna = Seq("AUGCGAUCGT", "Mixed RNA")
        # Это будет UNKNOWN из-за смешения T и U
        assert mixed_rna.get_alphabet_type() == "UNKNOWN"

    def test_protein_sequences(self):
        """Тест белковых последовательностей."""
        protein = Seq("MKAILVVLLYTRI", "Normal protein")
        assert protein.get_alphabet_type() == "PROTEIN"

        # Белок с стоп-кодоном
        protein_with_stop = Seq("MKAIL*", "Protein with stop")
        assert protein_with_stop.get_alphabet_type() == "PROTEIN"


    def test_reverse_complement(self):
        """Тест обратно-комплементарных последовательностей."""
        dna = Seq("ATGC", "Test DNA")
        rev_comp = dna.reverse_complement()
        assert rev_comp.sequence == "GCAT"

        # Тест с неоднозначными символами
        ambig_dna = Seq("ATGCRYSW", "Ambiguous DNA")
        rev_comp = ambig_dna.reverse_complement()
        expected = "WSYRGCAT"  # W↔S, Y↔R
        assert rev_comp.sequence == expected

    def test_translation(self):
        """Тест трансляции."""
        # Простая трансляция
        dna = Seq("ATGTAG", "Simple coding")
        protein = dna.translate()
        assert protein.sequence == "M*"  # ATG=M, TAG=*

        # Длинная последовательность
        long_dna = Seq("ATGAAACGCATTAGC", "Long coding")
        protein = long_dna.translate()
        assert len(protein.sequence) == 5  # 15 нуклеотидов = 5 аминокислот

    def test_composition_analysis(self):
        """Тест анализа состава."""
        seq = Seq("AATTGGCC", "Balanced")
        comp = seq.get_composition()
        assert comp == {'A': 2, 'T': 2, 'G': 2, 'C': 2}

        # Тест с белком
        protein = Seq("AAAAACCCCC", "Protein test")
        comp = protein.get_composition()
        assert comp == {'A': 5, 'C': 5}

    def test_formatting(self):
        """Тест форматирования вывода."""
        # Короткая последовательность
        short_seq = Seq("ATGC", "Short")
        formatted = str(short_seq)
        assert ">Short" in formatted
        assert "ATGC" in formatted

        # Длинная последовательность (тест переноса строк)
        long_seq = Seq("A" * 100, "Long sequence")
        formatted = str(long_seq)
        lines = formatted.split('
')
        assert len(lines) > 2  # Заголовок + минимум 2 строки последовательности

    def test_equality(self):
        """Тест сравнения объектов."""
        seq1 = Seq("ATGC", "Test1")
        seq2 = Seq("ATGC", "Test1")
        seq3 = Seq("ATGC", "Test2")  # Другой заголовок
        seq4 = Seq("TTTT", "Test1")  # Другая последовательность

        assert seq1 == seq2
        assert seq1 != seq3
        assert seq1 != seq4
        assert seq1 != "not a seq object"

    def test_error_handling(self):
        """Тест обработки ошибок."""
        # Пустая последовательность
        with pytest.raises(SequenceError):
            Seq("", "Empty")

        # Неправильный тип
        with pytest.raises(SequenceError):
            Seq(123, "Wrong type")


        # Обратная комплементарность для РНК
        rna = Seq("AUGC", "RNA")
        with pytest.raises(SequenceError):
            rna.reverse_complement()

    @pytest.mark.parametrize("sequence,expected_type", [
        ("ATGC", "DNA"),
        ("AUGC", "RNA"),
        ("MKAIL", "PROTEIN"),
        ("ATGCU", "UNKNOWN"),  # Смесь ДНК и РНК
        ("123", "UNKNOWN"),    # Недопустимые символы
    ])
    def test_sequence_type_detection(self, sequence, expected_type):
        """Параметризованный тест определения типа последовательности."""
        if expected_type == "UNKNOWN" and any(char in sequence for char in "123"):
            # Этот случай должен вызвать ошибку при создании
            with pytest.raises(SequenceError):
                Seq(sequence, "Test")
        else:
            seq = Seq(sequence, "Test")
            assert seq.get_alphabet_type() == expected_type

    def test_performance_large_sequences(self):
        """Тест производительности с большими последовательностями."""
        # Большая ДНК последовательность
        large_dna = Seq("A" * 10000 + "T" * 10000 + "G" * 10000 + "C" * 10000, "Large DNA")

        # Основные операции должны работать быстро
        assert len(large_dna) == 40000
        assert large_dna.get_alphabet_type() == "DNA"
        assert large_dna.get_gc_content() == 50.0

        # Состав должен быть корректным
        comp = large_dna.get_composition()
        assert comp == {'A': 10000, 'T': 10000, 'G': 10000, 'C': 10000}


if __name__ == '__main__':
    pytest.main([__file__])
