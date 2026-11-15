import Header from '../components/layout/Header';
import styles from './../globalStyles/text.module.css'

const Error = () => {
    return (
        <div>
            <Header />
            <p className={styles.title}>Page not found</p>
        </div>
    );
}

export default Error;