import { useEffect, useRef, useState } from 'react';
import styles from './css/styles.module.css'
import formStyles from './../globalStyles/form.module.css'
import wrapperStyles from './../globalStyles/wrappers.module.css'
import textStyles from './../globalStyles/text.module.css'
import PredTable from "../components/table/PredTable";
import BackButton from '../components/styledComponents/BackButton';
import ReagentSelector from '../components/styledComponents/ReagentSelector';

const ProcessMassValue = () => {
  const [comp, setComp] = useState();
  const [score, setScore] = useState();
  const [uComp, setUComp] = useState();
  const [PPMs, setPMMs] = useState();

  const [compDisplay, setCompDisplay] = useState();
  const [predDisplay, setPredDisplay] = useState();
  const [PPMDisplay, setPPMDisplay] = useState();

  const [simpleMode, setSimpleMode] = useState(false);
  const [cri1, setCr1] = useState(false);
  const [cri2, setCr2] = useState(false);
  const [cri3, setCr3] = useState(false);
  const [cri4, setCr4] = useState(false);

  const [criEnc, setCriEnc] = useState();

  const [showResults, setShowResults] = useState(false);
  const [loading, setLoading] = useState(false);
  const [errorMsg, setErrorMsg] = useState(null);
  const massRef = useRef(null);
  const dropdownRef = useRef(null);

  const submitHandler = (e) => {
    e.preventDefault();
    setErrorMsg(null);

    if (massRef.current.value.trim() == '' || massRef.current.value == null) {
      setErrorMsg("Enter an mz value");
      return;
    }

    setLoading(true);

    let hasErr = false;
    let statusCode;

    fetch('http://localhost:18080/predict', {
      method: 'POST',
      body: JSON.stringify({
        mz_array: [parseFloat(massRef.current.value) || -1],
        reagent_ion: dropdownRef.current.value || ''
      }),
      headers: {
        'Content-Type': 'application/json'
      }
    }).then((res) => {
      if (!res.ok) {
	statusCode = res.status;
	hasErr = true;
	return res.text();
      }

      return res.json();
    }).then((data) => {
      if (hasErr) {
	setErrorMsg(data);
	console.log(status + ": " + data);
	setLoading(false);
	return;
      }

      setLoading(false);
      
      data.compounds = data.compounds;
      data.uCompounds = data.uCompounds;
      data.ppms = data.ppms;
      data.scores = data.scores;

      setComp(data.compounds);
      setScore(data.scores);
      setUComp(data.uCompounds);
      setCriEnc(data.critereaEncodings);
      setCompDisplay(data.compounds);
      setPredDisplay(data.scores);
      setPPMDisplay(data.ppms)
      setPMMs(data.ppms);
      setShowResults(true);
      setSimpleMode(false);
    }).catch((err) => {
      setLoading(false);
      console.log(err);
      setErrorMsg("An error occured. Please try again. ");
    });
  }

  const showFormHandler = () => {
    setShowResults(false);
  }

  const doesPassCrit = (a, b) => {
    for (let i = 0; i < a.length; i++) {
      if (a[i] == 1 && b[i] == 0) {
        return false;
      }
    }

    return true;
  }

  useEffect(() => {
    let displayC = [];
    let displayP = score;
    let displayPPM = PPMs;

    if (!simpleMode) {
      displayC = comp;
    } else {
      displayC = uComp;
    }

    const enc = [Number(cri1), Number(cri2), Number(cri3), Number(cri4)];

    let isZero = true;
    for (const el of enc) {
      if (el !== 0) {
        isZero = false;
      }
    }

    if (!isZero) {
      let c = [];
      let p = [];
      let ppm = []

      for (let i = 0; i < criEnc.length; i++) {
        if (doesPassCrit(enc, criEnc[i])) {
          c.push(displayC[i]);
          p.push(displayP[i]);
          ppm.push(displayPPM[i]);
        }
      }
      
      displayC = c;
      displayP = p;
      displayPPM = ppm;
    }

    setCompDisplay(displayC);
    setPredDisplay(displayP);
    setPPMDisplay(displayPPM);
  }, [cri1, cri2, cri3, cri4, simpleMode, PPMs, comp, score, uComp, criEnc]);

  return (
    <div>
      <div className={wrapperStyles.three}>
        <div className={wrapperStyles['three-small']}>
          {showResults && <BackButton onClick={showFormHandler}>Back</BackButton>}
        </div>
        <div className={wrapperStyles['three-big']}>
          <p className={textStyles.title}>SoarAI 1.0</p>
        </div>
        <div className={wrapperStyles['three-small']}></div>
      </div>
      
      {!showResults &&
       <div className={wrapperStyles.center}>
	 <form onSubmit={submitHandler}>
	   <ReagentSelector ref={dropdownRef} />
	   <input type='type' placeholder='0.00' name='mass' ref={massRef}/>
	   <button type='submit'>Get Predictions</button>
	 </form>
       </div>
      }

      {errorMsg &&
       <div className={wrapperStyles.center}>
         <p className={formStyles.error}>{errorMsg}</p>    
       </div>}

      {loading &&
       <div className={wrapperStyles.center}>
         <p>Loading...</p>    
       </div>}

      {showResults &&
       <div className={styles['table-wrap']}>
	 <div className={styles['check-wrap']}>
	   <div className={styles.checkboxes}>
             <div>
               <input type='checkbox' checked={simpleMode} onChange={(e) => { setSimpleMode((prevState) => !prevState) }} />
               <label>Polyatomics</label>
             </div>
             
             <div>
               <input type='checkbox' checked={cri1} onChange={(e) => { setCr1((prevState) => !prevState) }} />
               <label>Even number of hydrogen in compounds with (NH4)</label>
             </div>

             <div>
               <input type='checkbox' checked={cri2} onChange={(e) => { setCr2((prevState) => !prevState) }} />
               <label>Hydrogen {"<="} ((2 * (# of Carbon)) + 2) in compounds with (NH4)</label>
             </div>

             <div>
               <input type='checkbox' checked={cri3} onChange={(e) => { setCr3((prevState) => !prevState) }} />
               <label>Odd number of Hydrogen in compounds with only Carbon, Hydrogen, and Oxygen</label>
             </div>
             
             <div>
               <input type='checkbox' checked={cri4} onChange={(e) => { setCr4((prevState) => !prevState) }} />
               <label>Hydrogen {"<="} ((2 * (# of Carbon)) + 3) in compounds with only Carbon, Hydrogen, and Oxygen</label>
             </div>
	   </div>
	 </div>

	 <PredTable compounds={compDisplay} preds={predDisplay} ppms={PPMDisplay} />
       </div>}

    </div>
  );
}

export default ProcessMassValue;
