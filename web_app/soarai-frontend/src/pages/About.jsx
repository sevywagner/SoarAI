import wrapperStyles from './../globalStyles/wrappers.module.css'
import styles from './css/styles.module.css'

const About = () => {
    return (
        <div>
            <div className={styles.thuggin}>
                <h1>He said that he on demon time... I said I'm on O(log(n))</h1>
                <img src={require('./../media/a.jpg')} height={500} />
            </div>
        </div>
    );
}

export default About;