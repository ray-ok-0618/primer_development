import streamlit as st

st.title("FASTA対応プライマー候補探索ツール（TM・GC・完全一致重視・閾値調整可能）")

iupac_dict = {
    frozenset(['A']): 'A',
    frozenset(['C']): 'C',
    frozenset(['G']): 'G',
    frozenset(['T']): 'T',
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['G', 'C']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C']): 'M',
    frozenset(['A', 'C', 'G']): 'V',
    frozenset(['A', 'C', 'T']): 'H',
    frozenset(['A', 'G', 'T']): 'D',
    frozenset(['C', 'G', 'T']): 'B',
    frozenset(['A', 'C', 'G', 'T']): 'N'
}

def get_iupac(bases):
    return iupac_dict.get(frozenset(bases), 'N')

def calc_base_frequencies(sequences, pos):
    bases = [seq[pos] for seq in sequences if seq[pos] in ['A', 'C', 'G', 'T']]
    total = len(bases)
    if total == 0:
        return {}
    freq = {}
    for b in ['A', 'C', 'G', 'T']:
        freq[b] = bases.count(b) / total
    return freq

def calc_tm_wallace(seq):
    seq = seq.upper()
    A = seq.count('A')
    T = seq.count('T')
    G = seq.count('G')
    C = seq.count('C')
    return 2 * (A + T) + 4 * (G + C)

def calc_gc_content(seq):
    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    length = len(seq.replace('-', ''))
    if length == 0:
        return 0
    return (gc_count / length) * 100

def read_fasta_alignment(lines):
    sequences = []
    seq_names = []
    current_seq = []

    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if current_seq:
                sequences.append(list(''.join(current_seq)))
                current_seq = []
            seq_names.append(line)
        else:
            current_seq.append(line)
    if current_seq:
        sequences.append(list(''.join(current_seq)))
    return seq_names, sequences

# ユーザーが閾値をスライダーで設定可能に
st.sidebar.header("プライマー条件設定")
min_tm = st.sidebar.slider("最小Tm値 (℃)", 0, 100, 55)
max_tm = st.sidebar.slider("最大Tm値 (℃)", 0, 100, 65)
min_gc = st.sidebar.slider("最小GC含有率 (%)", 0, 100, 40)
max_gc = st.sidebar.slider("最大GC含有率 (%)", 0, 100, 60)
min_len = st.sidebar.slider("最小塩基長", 0, 100, 20)
max_len = st.sidebar.slider("最大塩基長", 0, 100, 30)

def analyze_block(sequences, block_num=1):
    if len(sequences) == 0:
        st.write(f"--- ブロック{block_num}：配列なし ---")
        return

    length_set = set(len(seq) for seq in sequences)
    if len(length_set) != 1:
        st.warning(f"ブロック{block_num}: 配列長が揃っていません。処理を中止します。")
        return

    seq_len = length_set.pop()
    st.write(f"### ブロック {block_num}: 配列数={len(sequences)}, 配列長={seq_len}")

    max_rates = []
    consensus_seq = []

    for i in range(seq_len):
        freq = calc_base_frequencies(sequences, i)
        if not freq:
            consensus_seq.append('N')
            max_rates.append(0)
            continue
        max_freq = max(freq.values())
        max_rates.append(max_freq)
        if max_freq == 1.0:
            # 完全一致は単一塩基
            base = max(freq, key=freq.get)
            consensus_seq.append(base)
        else:
            # 完全一致でない場合は存在する塩基すべてで混合塩基
            bases_present = [b for b, r in freq.items() if r > 0]
            consensus_seq.append(get_iupac(bases_present))

    consensus_str = ''.join(consensus_seq)

    threshold = 0.8

    candidates = []
    max_fullmatch_count = 0

    for window_size in range(min_len, max_len + 1):
        for start in range(seq_len - window_size + 1):
            window_rates = max_rates[start:start + window_size]
            high_rate_count = sum(r >= threshold for r in window_rates)
            if high_rate_count >= window_size * 0.9:
                # ギャップ含む配列除外
                has_gap = any('-' in seq[start:start + window_size] for seq in sequences)
                if has_gap:
                    continue

                primer_seq = consensus_str[start:start + window_size]
                tm = calc_tm_wallace(primer_seq)
                gc = calc_gc_content(primer_seq)

                if min_tm <= tm <= max_tm and min_gc <= gc <= max_gc:
                    fullmatch_count = sum(1 for rate in window_rates if rate == 1.0)

                    if fullmatch_count > max_fullmatch_count:
                        max_fullmatch_count = fullmatch_count
                        candidates = [(start + 1, start + window_size, primer_seq, tm, gc, fullmatch_count)]
                    elif fullmatch_count == max_fullmatch_count:
                        candidates.append((start + 1, start + window_size, primer_seq, tm, gc, fullmatch_count))

    if candidates:
        st.subheader(f"プライマー候補領域 ブロック{block_num}（開始-終了 : 配列 (Tm℃, GC%, 完全一致塩基数)）")
        for start_pos, end_pos, seq, tm, gc, fullmatch_count in candidates:
            st.text(f"{start_pos}-{end_pos}: {seq} (Tm={tm:.1f}℃, GC={gc:.1f}%, 完全一致塩基数={fullmatch_count})")
    else:
        st.write(f"ブロック{block_num}で条件に合うプライマー候補領域が見つかりませんでした。")

uploaded_file = st.file_uploader("FASTA形式アライメントファイルをアップロード（utf-8）", type=["fa", "fasta", "txt"])

if uploaded_file is not None:
    lines = uploaded_file.read().decode('utf-8').splitlines()
    seq_names, sequences = read_fasta_alignment(lines)

    if len(sequences) == 0:
        st.error("配列が読み込めませんでした。FASTA形式ファイルを確認してください。")
    else:
        length_set = set(len(seq) for seq in sequences)
        if len(length_set) != 1:
            st.error("配列長が揃っていません。アラインメントされたFASTAファイルをアップロードしてください。")
        else:
            analyze_block(sequences, 1)
