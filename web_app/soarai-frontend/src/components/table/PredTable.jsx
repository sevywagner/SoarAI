import { useState } from 'react';
import styles from './css/table.module.css';
import TableItem from './TableItem';

const PredTable = ({ compounds, preds, ppms }) => {
    return (
        <>
            <div className={styles.center}>
                <table border='true' className={styles.table}>
                    <tr>
                        <th className={styles.big}>Prediction</th>
                        <th className={styles.small}>Confidence</th>
                        <th className={styles.smaller}>PPM</th>
                    </tr>
                    {compounds.map((c, idx) => <TableItem key={Math.random()} compound={c} pred={preds[idx].toFixed(4)} ppm={ppms[idx].toFixed(4)} />)}
                </table>
            </div>
        </>
    );
}

export default PredTable;
