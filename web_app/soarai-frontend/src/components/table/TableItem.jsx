import styles from './css/table.module.css';

const TableItem = ({ compound, pred, ppm }) => {
    return (
        <tr>
            <td className={styles.big}>{compound}</td>
            <td className={styles.small}>{pred}%</td>
            <td className={styles.smaller}>{ppm}</td>
        </tr>
    );
}

export default TableItem;